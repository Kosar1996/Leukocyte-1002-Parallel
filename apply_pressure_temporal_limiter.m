function [state, fluid, limited, reason] = ...
    apply_pressure_temporal_limiter(state, fluid, old, z, par)
    limited = false;
    reason = '';

    if ~isfield(par, 'pressureLimiterEnabled') || ~par.pressureLimiterEnabled
        return;
    end
    if isfield(par, 'useFull2DFluid') && par.useFull2DFluid
        reason = 'pressure limiter skipped for full-2D Stokes fluid output';
        return;
    end
    if ~isfield(old, 'p') || ~isfield(old, 'pPrev') || ~isfield(state, 'p')
        return;
    end

    pNew = state.p(:);
    pOld = old.p(:);
    pPrev = old.pPrev(:);
    if numel(pNew) ~= numel(pOld) || numel(pOld) ~= numel(pPrev)
        return;
    end
    if ~all(isfinite(pNew)) || ~all(isfinite(pOld)) || ~all(isfinite(pPrev))
        return;
    end

    N = numel(pNew);
    idx = 2:(N-1);
    if isempty(idx)
        idx = 1:N;
    end

    maxStepAbs = 25;
    if isfield(par, 'pressureLimiterMaxStepAbs') && ...
            isfinite(par.pressureLimiterMaxStepAbs) && par.pressureLimiterMaxStepAbs > 0
        maxStepAbs = par.pressureLimiterMaxStepAbs;
    end

    dpNew = pNew - pOld;
    if max(abs(pOld(idx))) < maxStepAbs && max(abs(pPrev(idx))) < maxStepAbs
        return;
    end

    dropRatio = max((-dpNew(idx)) / maxStepAbs);

    if dropRatio <= 1
        return;
    end

    alpha = min(1, 1 / dropRatio);
    pLimited = pOld + alpha * dpNew;
    pLimited(1) = par.pIn;
    pLimited(end) = par.pOut;

    state.p = pLimited;
    [Q, ~, ~, tauL, tauE, uzL, uzE] = ...
        local_flux_and_shear(z, pLimited, state.deltaL, state.deltaE, ...
                             state.UwL, state.UwE, par);

    state.tauE = tauE;
    state.tauL = tauL;
    state.uzE = uzE;
    state.uzL = uzL;

    if ~isempty(fluid)
        fluid.p = pLimited;
        fluid.Q = Q;
        fluid.tauL = tauL;
        fluid.tauE = tauE;
        fluid.uzL = uzL;
        fluid.uzE = uzE;
        fluid.gap = state.deltaE - state.deltaL;
    end

    limited = true;
    reason = sprintf('limited max pressure drop from %.3e Pa to %.3e Pa', ...
        max(-dpNew(idx)), max(-(pLimited(idx) - pOld(idx))));
end