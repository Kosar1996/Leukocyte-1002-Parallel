function vq = interp_curve_values(zNodes, vNodes, zq)
% INTERP_CURVE_VALUES
% Position-aware Q4 boundary profile interpolation with endpoint clamping.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 28-35 (Boundary Endpoint Clamping Safeguard):
%    [OLD]: vq = safe_interp1_same_or_resample(zu, vs, zq, 'interp_curve_values');
%    Evaluating query points (zq) outside the active node range [min(zu), max(zu)] 
%    with un-clamped linear extrapolation generated unphysical radial bulges 
%    and velocity spikes at domain ends (z = zMin and z = zMax).
%    FIXED: Clamped interpolated query values (vq) at boundary limits zu(1) 
%    and zu(end) to enforce flat physical profile continuity beyond boundary nodes.
% -------------------------------------------------------------------------

    zNodes = zNodes(:);
    vNodes = vNodes(:);
    zq = zq(:);

    [zs, idx] = sort(zNodes);
    vs = vNodes(idx);
    [zu, ~, ic] = unique(zs);
    if numel(zu) < numel(zs)
        vs = accumarray(ic, vs, [], @mean);
    end

    if numel(zu) == 1
        vq = vs(1) * ones(size(zq));
    else
        % [OLD]: vq = safe_interp1_same_or_resample(zu, vs, zq, 'interp_curve_values');
        vq = safe_interp1_same_or_resample(zu, vs, zq, 'interp_curve_values');
        
        % Boundary endpoint clamping safeguard
        vq(zq < zu(1))   = vs(1);
        vq(zq > zu(end)) = vs(end);
    end
end