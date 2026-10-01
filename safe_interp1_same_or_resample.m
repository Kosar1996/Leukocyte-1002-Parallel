function vq = safe_interp1_same_or_resample(x, v, xq, label)
%SAFE_INTERP1_SAME_OR_RESAMPLE Robust linear interpolation for coupled grids.
% MATLAB interp1 requires numel(x)==numel(v). During the staged 2D MAC
% post-solve, some vectors can live on the solid-interface nodes while x is
% the fluid coupling grid. If lengths differ, resample v over the same
% physical interval as x instead of crashing with "X and V must be of the
% same length." This keeps the accepted solid step usable and makes the
% mismatch explicit when debugVerbose is enabled upstream.

    if nargin < 4 || isempty(label)
        label = 'unnamed'; %#ok<NASGU>
    end

    x = x(:);
    v = v(:);
    xq = xq(:);

    if isempty(xq)
        vq = zeros(0,1);
        return;
    end
    if isempty(v)
        vq = nan(size(xq));
        return;
    end
    if isempty(x)
        vq = v(1) * ones(size(xq));
        return;
    end

    if numel(x) ~= numel(v)
        if numel(v) == 1
            vq = v(1) * ones(size(xq));
            return;
        end
        x = linspace(min(x), max(x), numel(v)).';
    end

    good = isfinite(x) & isfinite(v);
    x = x(good);
    v = v(good);

    if isempty(x)
        vq = nan(size(xq));
        return;
    end

    [xs, idx] = sort(x);
    vs = v(idx);
    [xu, ~, ic] = unique(xs);
    if numel(xu) < numel(xs)
        vs = accumarray(ic, vs, [], @mean);
    end

    if numel(xu) == 1
        vq = vs(1) * ones(size(xq));
    else
        vq = interp1(xu, vs, xq, 'linear', 'extrap');
    end
end