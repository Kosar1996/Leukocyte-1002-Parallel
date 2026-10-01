function hybrid = softlube_build_hybrid_gap1d_exterior2d_mesh(out, statePlot, fluidPlot, opts)
%SOFTLUBE_BUILD_HYBRID_GAP1D_EXTERIOR2D_MESH Build a hybrid mesh summary.
%
% The hybrid layout is:
%   - 1D lubrication mesh along the narrow gap interval.
%   - 2D exterior reservoir grids upstream and downstream of that interval.
%
% This helper constructs the hybrid mesh geometry for plotting/prototyping.
% When called with a runtime hybrid fluid struct, it also records the final
% pressure/flux coupling metadata from the mixed-dimensional solve.

    if nargin < 4 || isempty(opts)
        opts = struct();
    end
    opts = hybrid_defaults(out, opts);

    zBase = out.z(:);
    rLBase = statePlot.deltaL(:);
    rEBase = statePlot.deltaE(:);

    if numel(rLBase) ~= numel(zBase) || numel(rEBase) ~= numel(zBase)
        error('statePlot.deltaL/deltaE must have the same length as out.z.');
    end

    gapZ = sort(opts.gapZ(:)).';
    gapZ(1) = max(gapZ(1), opts.zMin);
    gapZ(2) = min(gapZ(2), opts.zMax);
    if gapZ(2) <= gapZ(1)
        error('Expected a non-empty gap interval.');
    end

    zGap = zBase(zBase >= gapZ(1) - opts.zTol & zBase <= gapZ(2) + opts.zTol);
    zGap = unique([gapZ(:); zGap(:)]);
    if numel(zGap) < 3
        zGap = linspace(gapZ(1), gapZ(2), opts.NzGapFallback).';
    end

    rLGap = interp_curve_values_local(zBase, rLBase, zGap);
    rEGap = interp_curve_values_local(zBase, rEBase, zGap);
    rMidGap = 0.5 * (rLGap + rEGap);

    pGap = nan(size(zGap));
    if isfield(fluidPlot, 'p') && numel(fluidPlot.p) == numel(zBase)
        pGap = interp_curve_values_local(zBase, fluidPlot.p(:), zGap);
    elseif isfield(out, 'pHist') && ~isempty(out.pHist)
        pGap = interp_curve_values_local(zBase, out.pHist(:, end), zGap);
    end

    zUp = exterior_z_nodes_graded(opts.zMin, gapZ(1), ...
        opts.exteriorDzFar, opts.exteriorDzNear, opts.exteriorNearLength, 'clusterEnd');
    zDown = exterior_z_nodes_graded(gapZ(2), opts.zMax, ...
        opts.exteriorDzFar, opts.exteriorDzNear, opts.exteriorNearLength, 'clusterStart');

    upstream = body_fitted_exterior_block(zBase, rLBase, rEBase, opts.NrExterior, zUp);
    downstream = body_fitted_exterior_block(zBase, rLBase, rEBase, opts.NrExterior, zDown);

    hybrid = struct();
    hybrid.type = 'gap1d_exterior2d';
    hybrid.gapZ = gapZ;
    hybrid.rDomain = [opts.rMin, opts.rOuter];
    hybrid.options = opts;
    hybrid.gap1D = struct( ...
        'z', zGap, ...
        'rL', rLGap, ...
        'rE', rEGap, ...
        'rMid', rMidGap, ...
        'h', rEGap - rLGap, ...
        'p', pGap);
    hybrid.exterior2D = struct( ...
        'upstream', upstream, ...
        'downstream', downstream);
    hybrid.coupling = struct( ...
        'z', gapZ, ...
        'rL', interp_curve_values_local(zBase, rLBase, gapZ(:)), ...
        'rE', interp_curve_values_local(zBase, rEBase, gapZ(:)), ...
        'condition', 'match pressure and per-radian flux between 2D exterior and 1D gap');
    hybrid.originalMACQuality = original_mac_quality(fluidPlot);
    hybrid.runtime = runtime_hybrid_metadata(fluidPlot);
    hybrid.note = ['Hybrid mesh: 1D lubrication gap plus 2D exterior ', ...
        'reservoir grids. If runtime.available is true, the pressure field ', ...
        'came from the mixed-dimensional solver path.'];
end

function runtime = runtime_hybrid_metadata(fluidPlot)
    runtime = struct('available', false);
    if ~isstruct(fluidPlot) || ~isfield(fluidPlot, 'meshType') || ...
            ~strcmpi(fluidPlot.meshType, 'hybrid_gap1d_exterior2d')
        return;
    end

    runtime.available = true;
    if isfield(fluidPlot, 'hybrid')
        runtime.gapZ = fluidPlot.hybrid.gapZ;
        runtime.useExterior2DPressureVector = fluidPlot.hybrid.useExterior2DPressureVector;
        runtime.note = fluidPlot.hybrid.note;
    end
    if isfield(fluidPlot, 'pReduced')
        runtime.p1D = fluidPlot.pReduced(:);
    end
    if isfield(fluidPlot, 'pHybrid')
        runtime.pHybrid = fluidPlot.pHybrid(:);
    end
    if isfield(fluidPlot, 'QReduced')
        runtime.Q1D = fluidPlot.QReduced(:);
    end
    if isfield(fluidPlot, 'Q')
        runtime.QHybrid = fluidPlot.Q(:);
    end
end

function opts = hybrid_defaults(out, opts)
    if ~isfield(opts, 'gapZ') || isempty(opts.gapZ)
        opts.gapZ = [-0.2e-6, 4.2e-6];
    end
    if ~isfield(opts, 'zMin') || isempty(opts.zMin)
        if isfield(out, 'par') && isfield(out.par, 'zMin')
            opts.zMin = out.par.zMin;
        else
            opts.zMin = min(out.z(:));
        end
    end
    if ~isfield(opts, 'zMax') || isempty(opts.zMax)
        if isfield(out, 'par') && isfield(out.par, 'zMax')
            opts.zMax = out.par.zMax;
        else
            opts.zMax = max(out.z(:));
        end
    end
    if ~isfield(opts, 'rOuter') || isempty(opts.rOuter)
        opts.rOuter = 15e-6;
        if isfield(out, 'par') && isfield(out.par, 'global1DOuterRadius')
            opts.rOuter = out.par.global1DOuterRadius;
        elseif isfield(out, 'global1D') && isfield(out.global1D, 'r')
            opts.rOuter = max(out.global1D.r(:));
        elseif isfield(out, 'RPHist') && ~isempty(out.RPHist)
            opts.rOuter = max(out.RPHist(:), [], 'omitnan');
        end
    end
    if ~isfield(opts, 'rMin') || isempty(opts.rMin)
        opts.rMin = 0;
    end
    if ~isfield(opts, 'NrExterior') || isempty(opts.NrExterior)
        opts.NrExterior = 33;
    end
    if ~isfield(opts, 'exteriorDzNear') || isempty(opts.exteriorDzNear)
        opts.exteriorDzNear = 0.1e-6;
    end
    if ~isfield(opts, 'exteriorDzFar') || isempty(opts.exteriorDzFar)
        opts.exteriorDzFar = 0.75e-6;
    end
    if ~isfield(opts, 'exteriorNearLength') || isempty(opts.exteriorNearLength)
        opts.exteriorNearLength = 1.0e-6;
    end
    if ~isfield(opts, 'NzGapFallback') || isempty(opts.NzGapFallback)
        opts.NzGapFallback = 81;
    end
    opts.zTol = max(1e-15, 100 * eps(max(abs([opts.zMin, opts.zMax, opts.gapZ(:).']))));
end

function block = body_fitted_exterior_block(zBase, rLBase, rEBase, Nr, zNodes)
    eta = exterior_eta_nodes(Nr);
    zNodes = zNodes(:);
    rL = interp_curve_values_local(zBase, rLBase, zNodes);
    rE = interp_curve_values_local(zBase, rEBase, zNodes);
    R = rL(:).' + eta(:) .* (rE(:).' - rL(:).');
    Z = repmat(zNodes(:).', numel(eta), 1);

    block = struct();
    block.r = [];
    block.z = zNodes(:);
    block.R = R;
    block.Z = Z;
    block.rInner = rL(:);
    block.rOuter = rE(:);
    block.eta = eta(:);
    block.solidNodeMask = false(size(R));
    block.fluidNodeMask = true(size(R));
end

function eta = exterior_eta_nodes(Nr)
    if Nr < 2
        eta = 0.5;
        return;
    end
    s = linspace(0, 1, Nr).';
    beta = 1.8;
    left = 0.5 * (2*s).^beta;
    right = 1 - 0.5 * (2*(1-s)).^beta;
    eta = zeros(size(s));
    eta(s <= 0.5) = left(s <= 0.5);
    eta(s > 0.5) = right(s > 0.5);
    eta(1) = 0;
    eta(end) = 1;
end

function z = exterior_z_nodes_graded(a, b, dzFar, dzNear, nearLength, clusterSide)
    if b < a
        z = zeros(0, 1);
        return;
    end
    if abs(b - a) < 1e-18
        z = [a; b];
        return;
    end

    nearLength = min(max(nearLength, 0), b - a);
    switch clusterSide
        case 'clusterEnd'
            zFar = uniform_segment(a, b - nearLength, dzFar);
            zNear = uniform_segment(b - nearLength, b, dzNear);
        case 'clusterStart'
            zNear = uniform_segment(a, a + nearLength, dzNear);
            zFar = uniform_segment(a + nearLength, b, dzFar);
        otherwise
            error('Unknown clusterSide: %s', clusterSide);
    end
    z = unique([a; b; zFar(:); zNear(:)]);
end

function z = uniform_segment(a, b, dz)
    if b <= a
        z = a;
        return;
    end
    n = max(1, ceil((b - a) / dz));
    z = linspace(a, b, n + 1).';
end

function q = original_mac_quality(fluidPlot)
    q = struct('available', false);
    if ~isfield(fluidPlot, 'meshF') || ~isfield(fluidPlot.meshF, 'Rp') || ...
            ~isfield(fluidPlot.meshF, 'Zp') || isempty(fluidPlot.meshF.Rp)
        return;
    end

    Rp = fluidPlot.meshF.Rp;
    Nr = size(Rp, 1);
    if Nr < 2
        return;
    end

    etaF = linspace(0, 1, Nr + 1).';
    etaC = 0.5 * (etaF(1:end-1) + etaF(2:end));
    hC = (Rp(end, :).' - Rp(1, :).') / (etaC(end) - etaC(1));
    ratio = abs(diff(hC)) ./ max(abs(hC(1:end-1)), 1e-30);

    q.available = true;
    q.minCenterGap = min(hC);
    q.maxCenterGap = max(hC);
    q.maxAdjacentGapChangeRatio = max(ratio);
end

function vq = interp_curve_values_local(zNodes, vNodes, zq)
    zNodes = zNodes(:);
    vNodes = vNodes(:);
    zq = zq(:);
    good = isfinite(zNodes) & isfinite(vNodes);
    zNodes = zNodes(good);
    vNodes = vNodes(good);

    if isempty(zNodes)
        vq = nan(size(zq));
        return;
    end

    [zs, idx] = sort(zNodes);
    vs = vNodes(idx);
    [zu, ~, ic] = unique(zs);
    if numel(zu) < numel(zs)
        vs = accumarray(ic, vs, [], @mean);
    end

    if isscalar(zu)
        vq = vs(1) * ones(size(zq));
    else
        vq = interp1(zu, vs, zq, 'linear', 'extrap');
    end
end
