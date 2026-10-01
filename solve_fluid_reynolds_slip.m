function [fluid, ok, stopReason] = solve_fluid_reynolds_slip(z, old, state, par)
% Solve the Reynolds-type fluid problem for the current coupled iterate.
%
% old   = previously accepted time-step state
% state = current inner coupling iterate at t^{n+1}

    N  = par.NzFluid;
    dzFace = diff(z(:));
    dz = par.dz;

    rl     = state.deltaL;
    re     = state.deltaE;
    rl_old = old.deltaL;
    re_old = old.deltaE;

    ok = true;
    stopReason = '';

    h = re - rl;

    if any(h <= par.minGap)
        [hmin, imin] = min(h);

        ok = false;
        stopReason = sprintf(['Gap fell below minGap during fluid solve. ', ...
                              'z = %.6e m, hmin = %.6e m, minGap = %.6e m'], ...
                              z(imin), hmin, par.minGap);

        fluid = struct();
        fluid.p    = state.p;
        fluid.Q    = nan(par.NzFluid-1,1);
        fluid.tauL = nan(par.NzFluid,1);
        fluid.tauE = nan(par.NzFluid,1);
        fluid.uzL  = nan(par.NzFluid,1);
        fluid.uzE  = nan(par.NzFluid,1);
        fluid.gap  = h;
        return;
    end

    % Area evolution source term
    A    = 0.5 * (re.^2     - rl.^2);
    Aold = 0.5 * (re_old.^2 - rl_old.^2);
    Ssrc = -par.SsrcFactor * (A - Aold) / par.dt;

    % fprintf('fluid: max|re-re_old| = %.3e, max|Ssrc| = %.3e\n', ...
    % max(abs(re-re_old)), max(abs(Ssrc)));

    % Pressure is algebraic, not a history variable. Start each fluid solve
    % from the boundary pressure to avoid recycling stale pressure plateaus.
    p      = linspace(par.pIn, par.pOut, N).';
    p(1)   = par.pIn;
    p(end) = par.pOut;

    if ~isfield(par, 'useDirectFluidSolve') || par.useDirectFluidSolve
        [Q, dQdpL, dQdpR, ~, ~, ~, ~] = ...
            local_flux_and_shear(z, p, rl, re, state.UwL, state.UwE, par);
        [R, J] = fluid_residual_jacobian(Q, p, dQdpL, dQdpR, Ssrc, dz, par);

        dp = -J \ R;
        if all(isfinite(dp))
            p = p + dp;
            p(1) = par.pIn;
            p(end) = par.pOut;

            [Q, ~, ~, tauL, tauE, uzL, uzE] = ...
                local_flux_and_shear(z, p, rl, re, state.UwL, state.UwE, par);
            Rfinal = fluid_residual_only(Q, p, Ssrc, dz, par);
            directTol = max(par.tolNewtonFluid, 1e-12);

            if norm(Rfinal, inf) < directTol
                fluid = struct('p', p, ...
                               'Q', Q, ...
                               'tauL', tauL, ...
                               'tauE', tauE, ...
                               'uzL', uzL, ...
                               'uzE', uzE, ...
                               'gap', h);
                return;
            end
        end

        if isfield(par, 'debugVerbose') && par.debugVerbose
            fprintf('   direct fluid solve residual too large; falling back to Newton\n');
        end
    end

    converged = false;

    for newt = 1:par.maxNewtonFluid
        [Q, dQdpL, dQdpR, tauL, tauE, uzL, uzE] = ...
            local_flux_and_shear(z, p, rl, re, state.UwL, state.UwE, par);

        [R, J] = fluid_residual_jacobian(Q, p, dQdpL, dQdpR, Ssrc, dz, par);

        resNorm = norm(R, inf);
        if resNorm < par.tolNewtonFluid
            converged = true;
            break;
        end

        dp = -J \ R;

        % Backtracking line search
        alpha    = 1.0;
        accepted = false;
        pTry     = p;

        for ls = 1:10
            pCand      = p + alpha * dp;
            pCand(1)   = par.pIn;
            pCand(end) = par.pOut;

            [Qtry, ~, ~, ~, ~, ~, ~] = ...
                local_flux_and_shear(z, pCand, rl, re, state.UwL, state.UwE, par);

            Rtry = fluid_residual_only(Qtry, pCand, Ssrc, dz, par);

            if norm(Rtry, inf) < resNorm
                pTry = pCand;
                accepted = true;
                break;
            end

            alpha = 0.5 * alpha;
        end

        if ~accepted
            ok = false;
            stopReason = sprintf('Fluid Newton line search failed at Newton iteration %d.', newt);

            fluid = struct();
            fluid.p    = p;
            fluid.Q    = nan(par.NzFluid-1,1);
            fluid.tauL = nan(par.NzFluid,1);
            fluid.tauE = nan(par.NzFluid,1);
            fluid.uzL  = nan(par.NzFluid,1);
            fluid.uzE  = nan(par.NzFluid,1);
            fluid.gap  = h;
            return;
        end

        p = pTry;
    end

    if ~converged
        % One last residual check after the final iterate
        [Q, ~, ~, tauL, tauE, uzL, uzE] = ...
            local_flux_and_shear(z, p, rl, re, state.UwL, state.UwE, par);
        Rfinal = fluid_residual_only(Q, p, Ssrc, dz, par);

        if norm(Rfinal, inf) >= par.tolNewtonFluid
            ok = false;
            stopReason = sprintf('Fluid Newton did not converge within %d iterations.', par.maxNewtonFluid);

            fluid = struct();
            fluid.p    = p;
            fluid.Q    = nan(par.NzFluid-1,1);
            fluid.tauL = nan(par.NzFluid,1);
            fluid.tauE = nan(par.NzFluid,1);
            fluid.uzL  = nan(par.NzFluid,1);
            fluid.uzE  = nan(par.NzFluid,1);
            fluid.gap  = h;
            return;
        end
    end

    [Q, ~, ~, tauL, tauE, uzL, uzE] = ...
        local_flux_and_shear(z, p, rl, re, state.UwL, state.UwE, par);

    fluid = struct('p', p, ...
                   'Q', Q, ...
                   'tauL', tauL, ...
                   'tauE', tauE, ...
                   'uzL', uzL, ...
                   'uzE', uzE, ...
                   'gap', h);
end

function [R, J] = fluid_residual_jacobian(Q, p, dQdpL, dQdpR, S, dz, par)
    N = par.NzFluid;
    R = fluid_residual_only(Q, p, S, dz, par);
    if isfield(par, 'zGrid') && numel(par.zGrid) == N
        dzControl = global_1d_control_lengths(par.zGrid(:));
    else
        dzControl = dz * ones(N,1);
    end

    nTrip = 2 + 4 * max(N - 2, 0);
    rows = zeros(nTrip,1);
    cols = zeros(nTrip,1);
    vals = zeros(nTrip,1);
    ptr = 1;

    rows(ptr) = 1;
    cols(ptr) = 1;
    vals(ptr) = 1;
    ptr = ptr + 1;

    for i = 2:N-1
        dxc = dzControl(i);
        rows(ptr:ptr+3) = i;
        cols(ptr:ptr+3) = [i; i+1; i-1; i];
        vals(ptr:ptr+3) = [dQdpL(i) / dxc; ...
                           dQdpR(i) / dxc; ...
                          -dQdpL(i-1) / dxc; ...
                          -dQdpR(i-1) / dxc];
        ptr = ptr + 4;
    end

    rows(ptr) = N;
    cols(ptr) = N;
    vals(ptr) = 1;

    J = sparse(rows, cols, vals, N, N);
end