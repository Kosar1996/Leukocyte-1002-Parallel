function [Rfaces, etaFaces] = build_radial_face_grid_bodyfitted(rInner, rOuter, Nr, par)
    rInner = rInner(:).';
    rOuter = rOuter(:).';
    nCol = numel(rInner);
    Rfaces = nan(Nr + 1, nCol);
    etaFaces = nan(Nr + 1, nCol);

    for j = 1:nCol
        a = rInner(j);
        b = rOuter(j);
        if ~(isfinite(a) && isfinite(b) && b > a)
            Rfaces(:,j) = linspace(a, b, Nr + 1).';
            etaFaces(:,j) = linspace(0, 1, Nr + 1).';
            continue;
        end

        rFaces = radial_faces_for_interval_bodyfitted(a, b, Nr, par);
        Rfaces(:,j) = rFaces(:);
        etaFaces(:,j) = (rFaces(:) - a) / max(b - a, eps);
    end
end

function rFaces = radial_faces_for_interval_bodyfitted(a, b, Nr, par)
    rFaces = linspace(a, b, Nr + 1).';
    if ~isfield(par, 'bodyFittedRadialFineWindow') || ...
            numel(par.bodyFittedRadialFineWindow) ~= 2
        return;
    end

    fineWindow = sort(par.bodyFittedRadialFineWindow(:)).';
    if any(~isfinite(fineWindow)) || fineWindow(2) <= fineWindow(1)
        return;
    end

    fineWeight = 4.0;
    if isfield(par, 'bodyFittedRadialFineWeight') && ...
            isfinite(par.bodyFittedRadialFineWeight) && ...
            par.bodyFittedRadialFineWeight > 1
        fineWeight = par.bodyFittedRadialFineWeight;
    end

    fineStart = max(a, fineWindow(1));
    fineEnd = min(b, fineWindow(2));
    tol = max(100 * eps(max(abs([a, b, fineWindow, 1]))), 1e-15);
    if fineEnd <= fineStart + tol
        return;
    end

    segStart = [];
    segEnd = [];
    segWeight = [];
    if fineStart > a + tol
        segStart(end+1) = a;
        segEnd(end+1) = fineStart;
        segWeight(end+1) = 1;
    end
    segStart(end+1) = fineStart;
    segEnd(end+1) = fineEnd;
    segWeight(end+1) = fineWeight;
    if b > fineEnd + tol
        segStart(end+1) = fineEnd;
        segEnd(end+1) = b;
        segWeight(end+1) = 1;
    end

    segLength = segEnd - segStart;
    good = segLength > tol;
    segStart = segStart(good);
    segEnd = segEnd(good);
    segWeight = segWeight(good);
    segLength = segLength(good);
    if isempty(segLength)
        return;
    end

    nCells = allocate_weighted_segment_cells(segLength, segWeight, Nr);
    rFaces = segStart(1);
    for s = 1:numel(nCells)
        localFaces = linspace(segStart(s), segEnd(s), nCells(s) + 1);
        rFaces = [rFaces(:); localFaces(2:end).'];
    end

    if numel(rFaces) ~= Nr + 1 || any(diff(rFaces) <= 0)
        rFaces = linspace(a, b, Nr + 1).';
    else
        rFaces(1) = a;
        rFaces(end) = b;
    end
end

function nCells = allocate_weighted_segment_cells(segLength, segWeight, Nr)
    score = segLength(:) .* segWeight(:);
    raw = Nr * score / sum(score);
    nCells = max(1, floor(raw));

    while sum(nCells) < Nr
        [~, idx] = max(raw - nCells);
        nCells(idx) = nCells(idx) + 1;
    end

    while sum(nCells) > Nr
        canReduce = find(nCells > 1);
        if isempty(canReduce)
            break;
        end
        [~, localIdx] = min(raw(canReduce) - nCells(canReduce));
        idx = canReduce(localIdx);
        nCells(idx) = nCells(idx) - 1;
    end
end