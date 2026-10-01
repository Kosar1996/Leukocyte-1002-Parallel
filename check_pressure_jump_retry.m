function check_pressure_jump_retry(old, fluidTrial, z, par, stepIndex)
% Reject a time step when the pressure changes too abruptly.
%
% This function checks both:
%   1. the axial pressure vector fluidTrial.p used by the legacy diagnostics;
%   2. the native 2D MAC pressure field fluidTrial.P, when a previous 2D
%      field is available in old.P2DField.
%
% It intentionally throws an error. The outer try/catch reduces dt and
% retries the same physical time step.

    if ~isfield(par,'enablePressureJumpRetry') || ~par.enablePressureJumpRetry
        return;
    end

    if nargin >= 5 && isfield(par,'pressureJumpIgnoreFirstStep') && ...
            par.pressureJumpIgnoreFirstStep && stepIndex == 0
        return;
    end

    % -------------------------------
    % Axial pressure-vector check p(z)
    % -------------------------------
    if isfield(fluidTrial,'p') && ~isempty(fluidTrial.p)
        pNew = fluidTrial.p(:);

        if isfield(old,'p2D') && numel(old.p2D) == numel(pNew)
            pOld = old.p2D(:);
        elseif isfield(old,'p') && numel(old.p) == numel(pNew)
            pOld = old.p(:);
        elseif isfield(old,'pReduced') && numel(old.pReduced) == numel(pNew)
            pOld = old.pReduced(:);
        else
            pOld = [];
        end

        if ~isempty(pOld)
            dp = pNew - pOld;
            if isfield(par,'pressureJumpUseMidpointOnly') && par.pressureJumpUseMidpointOnly
                N = numel(pNew);
                jmid = round(N/2);
                hw = 2;
                if isfield(par,'pressureJumpMidHalfWidth') && isfinite(par.pressureJumpMidHalfWidth)
                    hw = max(0, round(par.pressureJumpMidHalfWidth));
                end
                jj = max(1,jmid-hw):min(N,jmid+hw);
                vals = dp(jj);
                vals = vals(isfinite(vals));
                if isempty(vals)
                    error('Pressure jump too large: nonfinite axial pressure jump detected');
                end
                jumpMeasure = abs(mean(vals));
            else
                vals = abs(dp(isfinite(dp)));
                if isempty(vals)
                    error('Pressure jump too large: nonfinite axial pressure jump detected');
                end
                jumpMeasure = max(vals);
            end

            absTol = 50;
            if isfield(par,'maxPressureJumpAbs') && isfinite(par.maxPressureJumpAbs)
                absTol = par.maxPressureJumpAbs;
            end
            jumpTol = pressure_jump_tolerance(absTol, par, pOld, pNew);

            if isfield(par,'debugVerbose') && par.debugVerbose
                fprintf('   pressure-jump axial check: maxDp = %.3e Pa, limit = %.3e Pa\n', ...
                    jumpMeasure, jumpTol);
            end

            jumpSafety = get_pressure_jump_safety_factor(par);
            atDtMin = is_at_minimum_dt(par);
            acceptAtDtMin = ~isfield(par,'acceptPressureJumpAtDtMin') || ...
                par.acceptPressureJumpAtDtMin;

            if jumpMeasure > jumpSafety*jumpTol
                if atDtMin && acceptAtDtMin
                    warning(['Pressure jump exceeded limit at dtMin, so the ', ...
                             'step is accepted: axial maxDp = %.6e Pa, ', ...
                             'allowed = %.6e Pa, safetyFactor = %.3f'], ...
                             jumpMeasure, jumpTol, jumpSafety);
                else
                    error(['Pressure jump too large: axial maxDp = %.6e Pa, ', ...
                           'allowed = %.6e Pa, safetyFactor = %.3f'], ...
                           jumpMeasure, jumpTol, jumpSafety);
                end
            end
        end
    end

    % -------------------------------
    % Native full 2D pressure-field check P(r,z)
    % -------------------------------
    if isfield(fluidTrial,'P') && ~isempty(fluidTrial.P) && ...
            isfield(old,'P2DField') && ~isempty(old.P2DField) && ...
            isequal(size(fluidTrial.P), size(old.P2DField))

        Pnew = fluidTrial.P;
        Pold = old.P2DField;
        dP = Pnew - Pold;
        vals = abs(dP(isfinite(dP)));
        if isempty(vals)
            error('Pressure jump too large: nonfinite 2D pressure jump detected');
        end
        jump2D = max(vals);

        absTol2D = 75;
        if isfield(par,'maxPressureJump2DAbs') && isfinite(par.maxPressureJump2DAbs)
            absTol2D = par.maxPressureJump2DAbs;
        elseif isfield(par,'maxPressureJumpAbs') && isfinite(par.maxPressureJumpAbs)
            absTol2D = par.maxPressureJumpAbs;
        end
        jumpTol2D = pressure_jump_tolerance(absTol2D, par, Pold(:), Pnew(:));

        if isfield(par,'debugVerbose') && par.debugVerbose
            fprintf('   pressure-jump 2D check: maxDp2D = %.3e Pa, limit = %.3e Pa\n', ...
                jump2D, jumpTol2D);
        end

        jumpSafety = get_pressure_jump_safety_factor(par);
        atDtMin = is_at_minimum_dt(par);
        acceptAtDtMin = ~isfield(par,'acceptPressureJumpAtDtMin') || ...
            par.acceptPressureJumpAtDtMin;

        if jump2D > jumpSafety*jumpTol2D
            if atDtMin && acceptAtDtMin
                warning(['2D pressure jump exceeded limit at dtMin, so the ', ...
                         'step is accepted: maxDp2D = %.6e Pa, ', ...
                         'allowed = %.6e Pa, safetyFactor = %.3f'], ...
                         jump2D, jumpTol2D, jumpSafety);
            else
                error(['Pressure jump too large: 2D maxDp = %.6e Pa, ', ...
                       'allowed = %.6e Pa, safetyFactor = %.3f'], ...
                       jump2D, jumpTol2D, jumpSafety);
            end
        end
    end
end

function jumpTol = pressure_jump_tolerance(absTol, par, pOld, pNew)
% Return a pressure-jump tolerance. If par.maxPressureJumpRel is finite,
% use it as an additional relative tolerance with absTol as the lower bound.
% If it is Inf or missing, use the absolute tolerance only.
    relTol = inf;
    if isfield(par,'maxPressureJumpRel') && isfinite(par.maxPressureJumpRel)
        relTol = par.maxPressureJumpRel;
    end

    if isfinite(relTol)
        pVals = [abs(pOld(isfinite(pOld))); abs(pNew(isfinite(pNew)))];
        if isempty(pVals)
            pScale = 1;
        else
            pScale = max(max(pVals), 1);
        end
        jumpTol = max(absTol, relTol*pScale);
    else
        jumpTol = absTol;
    end
end

function jumpSafety = get_pressure_jump_safety_factor(par)
% Multiplicative buffer for pressure-jump rejection. A value slightly above
% one prevents stopping for roundoff-level overshoots, e.g. 25.008 Pa when
% the limit is 25 Pa.
    jumpSafety = 1.05;
    if isfield(par,'pressureJumpSafetyFactor') && ...
            isfinite(par.pressureJumpSafetyFactor) && ...
            par.pressureJumpSafetyFactor >= 1
        jumpSafety = par.pressureJumpSafetyFactor;
    end
end

function tf = is_at_minimum_dt(par)
% True when the currently attempted step par.dt has already reached dtMin.
    tf = false;
    if isfield(par,'dt') && isfield(par,'dtMin') && ...
            isfinite(par.dt) && isfinite(par.dtMin) && par.dtMin > 0
        tf = par.dt <= par.dtMin*(1 + 1e-10);
    end
end