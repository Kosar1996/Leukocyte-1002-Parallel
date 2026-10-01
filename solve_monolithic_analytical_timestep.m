function stateNew = solve_monolithic_analytical_timestep(old, meshE, interfaceE, baseE, z, par)
% Fully coupled Newton step for the endothelium displacement and pressure.
%
% The accepted iterate must reduce the combined solid/fluid residual, and
% convergence requires both residuals to be below their absolute tolerances.
% minGap is only an admissibility constraint; the wall is never forced toward
% a prescribed gap.

    if isfield(par, 'useFsolveMono') && par.useFsolveMono
        stateNew = solve_monolithic_fsolve_timestep(old, meshE, interfaceE, baseE, z, par);
        return;
    end

    ndof = size(meshE.nodes,1)*2;
    N = numel(z);

    [fixDofs, fixVals] = solid_support_conditions(baseE, par.supportE);
    free = setdiff((1:ndof).', unique(fixDofs(:)));

    [Hr, Hz] = build_interface_interp_matrices(meshE, interfaceE, z); %#ok<ASGLU>

    u = old.uE;
    p = old.p;

    u(fixDofs) = fixVals;
    p(1) = par.pIn;
    p(end) = par.pOut;

    trustU = par.trustU0;
    [solidTarget, fluidTarget] = monolithic_residual_targets(par);

    bestMerit = inf;
    bestU = u;
    bestP = p;
    bestSolid = inf;
    bestFluid = inf;

    for newt = 1:par.maxNewtonMono

        [Rsolid, Rfluid, Kuu, Kup, Kpu, Kpp] = ...
            monolithic_analytical_residual_jacobian( ...
            u, p, old, meshE, interfaceE, z, par, ...
            free, fixDofs, fixVals, Hr, Hz);

        if newt == 1 && isfield(par, 'debugCoupledJacobian') && par.debugCoupledJacobian
            check_coupled_jacobian(u, p, old, meshE, interfaceE, z, par, ...
                free, fixDofs, fixVals, Kuu, Kup, Kpu, Kpp, Rsolid, Rfluid);
            error('Stopped after coupled Jacobian diagnostic.');
        end

        solidNorm = norm(Rsolid, inf);
        fluidNorm = norm(Rfluid, inf);
        merit = coupled_merit(solidNorm, fluidNorm, solidTarget, fluidTarget);

        [deltaE, ~] = monolithic_interface_kinematics( ...
            meshE, u, old.uE, interfaceE, z, par, Hr, Hz);
        gapMin = min(deltaE - old.deltaL);

        fprintf(['   Coupled Newton %d: solid=%.3e, fluid=%.3e, ', ...
                 'merit=%.3e, gap=%.3e, maxP=%.3e\n'], ...
                newt, solidNorm, fluidNorm, merit, gapMin, max(p));

        if merit < bestMerit
            bestMerit = merit;
            bestU = u;
            bestP = p;
            bestSolid = solidNorm;
            bestFluid = fluidNorm;
        end

        if solidNorm < solidTarget && fluidNorm < fluidTarget
            fprintf('   Coupled solid-fluid equilibrium converged\n');
            stateNew = build_coupled_state_from_up(u, p, old, meshE, interfaceE, z, par);
            return;
        end

        solidScale = max(solidNorm, solidTarget);
        fluidScale = max(fluidNorm, fluidTarget);

        Ju = par.uScaleMono;
        Jp = par.pScaleMono;

        J = [(Kuu * Ju) / solidScale, (Kup * Jp) / solidScale;
             (Kpu * Ju) / fluidScale, (Kpp * Jp) / fluidScale];

        R = [Rsolid / solidScale;
             Rfluid / fluidScale];

        dy = -J \ R;

        duFree = dy(1:numel(free)) * Ju;
        dpInt  = dy(numel(free)+1:end) * Jp;

        if ~all(isfinite(duFree)) || ~all(isfinite(dpInt))
            error('Coupled Newton produced non-finite increments.');
        end

        duMax = max(abs(duFree));
        if duMax > trustU
            stepScale = trustU / duMax;
        else
            stepScale = 1.0;
        end

        fprintf('      Newton increment: duMax=%.3e, dpMax=%.3e, trustU=%.3e\n', ...
                duMax, max(abs(dpInt)), trustU);

        accepted = false;
        alpha = 1.0;

        for ls = 1:par.lineSearchMaxMono

            uTrial = u;
            pTrial = p;

            uTrial(free) = uTrial(free) + alpha * stepScale * duFree;
            uTrial(fixDofs) = fixVals;

            pTrial(2:end-1) = pTrial(2:end-1) + alpha * stepScale * dpInt;
            pTrial(1) = par.pIn;
            pTrial(end) = par.pOut;

            try
                [deltaETrial, ~] = monolithic_interface_kinematics( ...
                    meshE, uTrial, old.uE, interfaceE, z, par, Hr, Hz);
                if any(deltaETrial - old.deltaL <= par.minGap)
                    alpha = 0.5 * alpha;
                    continue;
                end
                assert_solid_geometry_ok( ...
                    solid_geometry_quality(meshE, uTrial, 'endothelium'), par);

                [RsTrial, RfTrial] = monolithic_analytical_residual_only( ...
                    uTrial, pTrial, old, meshE, interfaceE, z, par, ...
                    free, fixDofs, fixVals);

                solidTrial = norm(RsTrial, inf);
                fluidTrial = norm(RfTrial, inf);
                meritTrial = coupled_merit(solidTrial, fluidTrial, solidTarget, fluidTarget);

                fprintf(['      coupled LS %d: alpha=%.3e, step=%.3e, ', ...
                         'solid=%.3e, fluid=%.3e, merit=%.3e\n'], ...
                        ls, alpha, stepScale, solidTrial, fluidTrial, meritTrial);

                if meritTrial < bestMerit
                    bestMerit = meritTrial;
                    bestU = uTrial;
                    bestP = pTrial;
                    bestSolid = solidTrial;
                    bestFluid = fluidTrial;
                end

                if meritTrial < merit
                    u = uTrial;
                    p = pTrial;
                    accepted = true;

                    if meritTrial < 0.25 * merit
                        trustU = min(2.0 * trustU, par.trustUMax);
                    elseif alpha < 0.25 || stepScale < 0.25
                        trustU = max(0.5 * trustU, par.trustUMin);
                    end

                    break;
                end

            catch ME
	                if ~(contains(ME.message, 'Negative or zero J') || ...
	                     contains(ME.message, 'Non-positive radius') || ...
	                     contains(ME.message, 'Element inverted') || ...
	                     contains(ME.message, 'Solid geometry guard failed') || ...
	                     contains(ME.message, 'Gap violates minGap'))
	                    rethrow(ME);
	                end
            end

            alpha = 0.5 * alpha;
        end

        if ~accepted
            trustU = 0.5 * trustU;
            if trustU <= par.trustUMin
                error(['Coupled Newton line search failed before equilibrium. ', ...
                       'best solid = %.3e, best fluid = %.3e, best merit = %.3e.'], ...
                       bestSolid, bestFluid, bestMerit);
            end
        end
    end

    error(['Coupled solid-fluid equilibrium did not converge within %d iterations. ', ...
           'best solid = %.3e, best fluid = %.3e, best merit = %.3e.'], ...
          par.maxNewtonMono, bestSolid, bestFluid, bestMerit);
end

function stateNew = solve_monolithic_fsolve_timestep(old, meshE, interfaceE, baseE, z, par)

    ndof = size(meshE.nodes,1)*2;

    [fixDofs, fixVals] = solid_support_conditions(baseE, par.supportE);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    [Hr, Hz] = build_interface_interp_matrices(meshE, interfaceE, z);

    u0 = old.uE;
    p0 = old.p;
    if ~all(isfinite(p0)) || numel(p0) ~= numel(z)
        p0 = linspace(par.pIn, par.pOut, numel(z)).';
    end

    if isfield(par, 'useMonoPredictor') && par.useMonoPredictor && ...
            isfield(old, 'uEPrev') && isfield(old, 'pPrev') && isfield(old, 'dtPrev')
        predScale = min(1.0, par.dt / max(old.dtPrev, eps));
        uPred = old.uE + predScale * (old.uE - old.uEPrev);
        pPred = old.p  + predScale * (old.p  - old.pPrev);
        pPred(1) = par.pIn;
        pPred(end) = par.pOut;

        [deltaPred, ~] = monolithic_interface_kinematics( ...
            meshE, uPred, old.uE, interfaceE, z, par, Hr, Hz);
        predGeometryOk = false;
        if all(isfinite(uPred))
            try
                assert_solid_geometry_ok( ...
                    solid_geometry_quality(meshE, uPred, 'endothelium predictor'), par);
                predGeometryOk = true;
            catch
                predGeometryOk = false;
            end
        end
        if all(isfinite(uPred)) && all(isfinite(pPred)) && ...
                all(deltaPred - old.deltaL > par.minGap) && predGeometryOk
            u0 = uPred;
            p0 = pPred;
        end
    end

    u0(fixDofs) = fixVals;
    p0(1) = par.pIn;
    p0(end) = par.pOut;

    Ju = par.uScaleMono;
    Jp = par.pScaleMono;

    y0 = [u0(free) / Ju;
          p0(2:end-1) / Jp];

    [solidTarget, fluidTarget] = monolithic_residual_targets(par);

    fun = @(y) monolithic_scaled_residual_from_y( ...
        y, old, meshE, interfaceE, z, par, free, fixDofs, fixVals, Hr, Hz, ...
        Ju, Jp, solidTarget, fluidTarget);

    algorithms = {'trust-region-dogleg'};
    if isfield(par, 'fsolveAlgorithms')
        algorithms = par.fsolveAlgorithms;
    end

    bestY = y0;
    bestSolid = inf;
    bestFluid = inf;
    bestExitflag = NaN;
    bestOutput = struct('iterations', 0);
    accepted = false;

    for attempt = 1:numel(algorithms)
        fsolveDisplay = 'iter';
        if isfield(par, 'fsolveDisplay')
            fsolveDisplay = par.fsolveDisplay;
        end
        printFsolve = ~strcmpi(fsolveDisplay, 'off') || ...
            (isfield(par, 'debugVerbose') && par.debugVerbose);

        opts = optimoptions('fsolve', ...
            'Algorithm', algorithms{attempt}, ...
            'SpecifyObjectiveGradient', true, ...
            'Display', fsolveDisplay, ...
            'MaxIterations', par.maxNewtonMono, ...
            'MaxFunctionEvaluations', 20 * numel(y0), ...
            'FunctionTolerance', 1e-10, ...
            'StepTolerance', 1e-14, ...
            'OptimalityTolerance', 1e-10);

        if isfield(par, 'useFsolveOutputFcn') && par.useFsolveOutputFcn
            opts = optimoptions(opts, 'OutputFcn', ...
                @(y, optimValues, state) stop_fsolve_when_coupled_converged( ...
                y, optimValues, state, old, meshE, interfaceE, z, par, ...
                free, fixDofs, fixVals, Ju, Jp, solidTarget, fluidTarget));
        end

        if printFsolve
            fprintf('   fsolve attempt %d using %s\n', attempt, algorithms{attempt});
        end
        [yTry, ~, exitflagTry, outputTry] = fsolve(fun, bestY, opts);

        [solidTry, fluidTry, gapTry] = coupled_residual_norms_from_y( ...
            yTry, old, meshE, interfaceE, z, par, free, fixDofs, fixVals, Ju, Jp);

        if printFsolve
            fprintf(['   fsolve attempt %d result: exitflag=%d, iterations=%d, ', ...
                     'solid=%.3e, fluid=%.3e, gap=%.3e\n'], ...
                    attempt, exitflagTry, outputTry.iterations, solidTry, fluidTry, gapTry);
        end

        if max(solidTry / solidTarget, fluidTry / fluidTarget) < ...
                max(bestSolid / solidTarget, bestFluid / fluidTarget)
            bestY = yTry;
            bestSolid = solidTry;
            bestFluid = fluidTry;
            bestExitflag = exitflagTry;
            bestOutput = outputTry;
        end

        if solidTry < solidTarget && fluidTry < fluidTarget && gapTry > par.minGap
            accepted = true;
            break;
        end
    end

    y = bestY;
    exitflag = bestExitflag;
    output = bestOutput;

    u = u0;
    p = p0;
    u(free) = y(1:numel(free)) * Ju;
    u(fixDofs) = fixVals;
    p(2:end-1) = y(numel(free)+1:end) * Jp;
    p(1) = par.pIn;
    p(end) = par.pOut;

    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshE, u, 'endothelium fsolve solution'), par);

    [Rsolid, Rfluid] = monolithic_analytical_residual_only( ...
        u, p, old, meshE, interfaceE, z, par, free, fixDofs, fixVals);

    solidNorm = norm(Rsolid, inf);
    fluidNorm = norm(Rfluid, inf);
    [deltaE, ~] = monolithic_interface_kinematics_value_only( ...
        meshE, u, old.uE, interfaceE, z, par);
    gapMin = min(deltaE - old.deltaL);

    if ~isfield(par, 'fsolveDisplay') || ~strcmpi(par.fsolveDisplay, 'off') || ...
            (isfield(par, 'debugVerbose') && par.debugVerbose)
        fprintf(['   fsolve exitflag=%d, iterations=%d, solid=%.3e, ', ...
                 'fluid=%.3e, gap=%.3e, maxP=%.3e\n'], ...
                exitflag, output.iterations, solidNorm, fluidNorm, gapMin, max(p));
    end

    if any(deltaE - old.deltaL <= par.minGap)
        error('Coupled fsolve solution violates minGap.');
    end

    if solidNorm >= solidTarget || fluidNorm >= fluidTarget
        error(['Coupled fsolve did not reach equilibrium. ', ...
               'solid = %.3e, fluid = %.3e.'], solidNorm, fluidNorm);
    end

    if ~accepted
        fprintf('   using best fsolve retry result before final residual check\n');
    end

    stateNew = build_coupled_state_from_up(u, p, old, meshE, interfaceE, z, par);
end


function [R, J] = monolithic_scaled_residual_from_y( ...
    y, old, meshE, interfaceE, z, par, free, fixDofs, fixVals, Hr, Hz, ...
    Ju, Jp, solidTarget, fluidTarget)

    ndof = size(meshE.nodes,1)*2;
    N = numel(z);

    u = old.uE;
    p = old.p;

    u(free) = y(1:numel(free)) * Ju;
    u(fixDofs) = fixVals;

    p(2:end-1) = y(numel(free)+1:end) * Jp;
    p(1) = par.pIn;
    p(end) = par.pOut;

    try
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshE, u, 'endothelium fsolve iterate'), par);
        if nargout > 1
            [Rsolid, Rfluid, Kuu, Kup, Kpu, Kpp] = ...
                monolithic_analytical_residual_jacobian( ...
                u, p, old, meshE, interfaceE, z, par, ...
                free, fixDofs, fixVals, Hr, Hz);

            J = [(Kuu * Ju) / solidTarget, (Kup * Jp) / solidTarget;
                 (Kpu * Ju) / fluidTarget, (Kpp * Jp) / fluidTarget];
        else
            [Rsolid, Rfluid] = monolithic_analytical_residual_only( ...
                u, p, old, meshE, interfaceE, z, par, ...
                free, fixDofs, fixVals);
        end

        R = [Rsolid / solidTarget;
             Rfluid / fluidTarget];

    catch ME
        if contains(ME.message, 'Negative or zero J') || ...
           contains(ME.message, 'Non-positive radius') || ...
           contains(ME.message, 'Element inverted') || ...
           contains(ME.message, 'Solid geometry guard failed') || ...
           contains(ME.message, 'Gap violates minGap')
            R = 1e12 * ones(numel(free) + N - 2, 1);
            if nargout > 1
                J = speye(numel(R), numel(y));
            end
        else
            rethrow(ME);
        end
    end
end

function [Rsolid, Rfluid, Kuu, Kup, Kpu, Kpp] = ...
    monolithic_analytical_residual_jacobian( ...
    u, p, old, meshE, interfaceE, z, par, ...
    free, fixDofs, fixVals, Hr, Hz)

    ndof = size(meshE.nodes,1)*2;
    N = numel(z);

    u(fixDofs) = fixVals;
    p(1) = par.pIn;
    p(end) = par.pOut;

    % Endothelium interface geometry used by the fluid solve.
    % If par.useExactDeformedInterfaceInMonolithic is true, deltaE and UwE
    % are evaluated at fixed fluid-grid positions using the deformed
    % interface z = Z + u_z, and HrUse/HUUse include the corresponding
    % geometric sensitivities.
    [deltaE, UwE, HrUse, HUUse] = monolithic_interface_kinematics( ...
        meshE, u, old.uE, interfaceE, z, par, Hr, Hz);

    deltaL = old.deltaL;
    h = deltaE - deltaL;

    if any(h <= par.minGap)
        error('Gap violates minGap.');
    end

    UwL = par.UwL * ones(N,1);

    sens = local_flux_shear_sensitivities(z, p, deltaL, deltaE, UwL, UwE, par);

    Q = sens.Q;
    tauE = sens.tauE;

    if isfield(par, 'debugVerbose') && par.debugVerbose
        fprintf('max |Q| = %.3e, max |diff(Q)/dz| = %.3e\n', ...
            max(abs(Q)), max(abs(diff(Q)/par.dz)));
    end

    trE.normal = p;
    trE.tangent = -tauE;
    if use_exact_interface_in_monolithic(par)
        trE.z = z;
    end

    Fext = zeros(ndof,1);

    [Fext, Kext, Bnormal, Btangent] = ...
        apply_interface_traction_sensitivity(meshE, u, Fext, interfaceE, trE);

    [Fint, Ktan] = assemble_finite_def_axisym(meshE, u, par);
    [Fvisc, Kvisc] = assemble_axisym_kelvin_voigt_viscous(meshE, u, old.uE, par);
    Fint = Fint + Fvisc;
    Ktan = Ktan + Kvisc;

    RsolidFull = Fint - Fext;
    Rsolid = RsolidFull(free);

    if isfield(par, 'debugVerbose') && par.debugVerbose
        fprintf('max Fint free = %.3e\n', norm(Fint(free),inf));
        fprintf('max Fext free = %.3e\n', norm(Fext(free),inf));
        fprintf('max Rsolid free = %.3e\n', norm(RsolidFull(free),inf));
        fprintf('max Rsolid full = %.3e\n', norm(RsolidFull,inf));
        fprintf('max Rsolid fixed = %.3e\n', norm(RsolidFull(fixDofs),inf));
    end

    re    = deltaE(:);
    rl    = deltaL(:);
    reOld = old.deltaE(:);
    rlOld = old.deltaL(:);

    A = 0.5*(re.^2 - rl.^2);
    Aold = 0.5*(reOld.^2 - rlOld.^2);

    Ssrc = -par.SsrcFactor*(A - Aold) / par.dt;

    Rfluid = zeros(N-2,1);

    for i = 2:N-1
        Rfluid(i-1) = (Q(i) - Q(i-1))/par.dz - Ssrc(i);
    end

    % -------------------------
    % Kpp = dRfluid/dp
    % -------------------------
    KppFull = sparse(N,N);

    for i = 2:N-1
        KppFull(i,i)   = KppFull(i,i)   + sens.dQdpL(i)   / par.dz;
        KppFull(i,i+1) = KppFull(i,i+1) + sens.dQdpR(i)   / par.dz;

        KppFull(i,i-1) = KppFull(i,i-1) - sens.dQdpL(i-1) / par.dz;
        KppFull(i,i)   = KppFull(i,i)   - sens.dQdpR(i-1) / par.dz;
    end

    Kpp = KppFull(2:end-1, 2:end-1);

    % -------------------------
    % Kpu = dRfluid/du
    % -------------------------
    Kp_delta = sparse(N-2,N);
    Kp_UwE   = sparse(N-2,N);

    for i = 2:N-1

        row = i-1;

        Kp_delta(row,i)   = Kp_delta(row,i)   + 0.5*sens.dQdb(i)   / par.dz;
        Kp_delta(row,i+1) = Kp_delta(row,i+1) + 0.5*sens.dQdb(i)   / par.dz;

        Kp_UwE(row,i)     = Kp_UwE(row,i)     + 0.5*sens.dQdUwE(i) / par.dz;
        Kp_UwE(row,i+1)   = Kp_UwE(row,i+1)   + 0.5*sens.dQdUwE(i) / par.dz;

        Kp_delta(row,i-1) = Kp_delta(row,i-1) - 0.5*sens.dQdb(i-1)   / par.dz;
        Kp_delta(row,i)   = Kp_delta(row,i)   - 0.5*sens.dQdb(i-1)   / par.dz;

        Kp_UwE(row,i-1)   = Kp_UwE(row,i-1)   - 0.5*sens.dQdUwE(i-1) / par.dz;
        Kp_UwE(row,i)     = Kp_UwE(row,i)     - 0.5*sens.dQdUwE(i-1) / par.dz;

        % source term contribution:
        % R = flux divergence - Ssrc
        % Ssrc = -(A-Aold)/dt
        % dR/dre = +re/dt
        Kp_delta(row,i) = Kp_delta(row,i) + par.SsrcFactor * re(i)/par.dt;
    end

    KpuFull = Kp_delta * HrUse + Kp_UwE * HUUse;
    Kpu = KpuFull(:,free);


        % % test numerical
        % % --- temporary numerical check for Kpu ---
        % % --- cheap numerical check for selected Kpu columns ---
        % epsU = 1e-12;
        % 
        % interfaceDofs = sort([2*interfaceE(:)-1; 2*interfaceE(:)]);
        % [~,loc] = ismember(interfaceDofs, free);
        % testCols = loc(loc > 0);
        % testCols = testCols(round(linspace(1,numel(testCols),min(20,numel(testCols)))));
        % 
        % [~, Rfluid0] = monolithic_analytical_residual_only( ...
        %     u, p, old, meshE, interfaceE, z, par, ...
        %     free, fixDofs, fixVals);
        % 
        % for jj = 1:numel(testCols)
        %     j = testCols(jj);
        % 
        %     uPert = u;
        %     uPert(free(j)) = uPert(free(j)) + epsU;
        % 
        %     [~, RfluidPert] = monolithic_analytical_residual_only( ...
        %         uPert, p, old, meshE, interfaceE, z, par, ...
        %         free, fixDofs, fixVals);
        % 
        %     Kpu_col_num = (RfluidPert - Rfluid0) / epsU;
        %     Kpu_col_ana = Kpu(:,j);
        % 
        %     err = norm(Kpu_col_ana - Kpu_col_num, inf) / ...
        %           max(norm(Kpu_col_num, inf), 1);
        % 
        %     fprintf('Kpu col %d rel error = %.3e\n', j, err);
        % end
        % % test end
        
    % -------------------------
    % Kuu and Kup
    % -------------------------
    Ip = speye(N);

    Ktau_p = -sens.dtauEdp;
    Ktau_u = -(sens.dtauEdb * HrUse + sens.dtauEdUwE * HUUse);

    KuuFull = Ktan - Kext - Btangent*Ktau_u;
    Kuu = KuuFull(free,free);

    KupFull = -(Bnormal*Ip + Btangent*Ktau_p);
    Kup = KupFull(free,2:end-1);

    % % ============================================================
    % % DEBUG: selected-column finite-difference check for Kuu
    % % ============================================================
    % epsU = 1e-12;
    % 
    % % test mostly interface DOFs, because these affect traction/gap most strongly
    % interfaceDofs = sort([2*interfaceE(:)-1; 2*interfaceE(:)]);
    % [~,loc] = ismember(interfaceDofs, free);
    % testCols = loc(loc > 0);
    % 
    % if isempty(testCols)
    %     testCols = round(linspace(1, numel(free), min(20,numel(free))));
    % else
    %     testCols = testCols(round(linspace(1,numel(testCols),min(20,numel(testCols)))));
    % end
    % 
    % [Rsolid0, ~] = monolithic_analytical_residual_only( ...
    %     u, p, old, meshE, interfaceE, z, par, ...
    %     free, fixDofs, fixVals);
    % 
    % for jj = 1:numel(testCols)
    %     j = testCols(jj);
    % 
    %     uPert = u;
    %     uPert(free(j)) = uPert(free(j)) + epsU;
    % 
    %     [RsolidPert, ~] = monolithic_analytical_residual_only( ...
    %         uPert, p, old, meshE, interfaceE, z, par, ...
    %         free, fixDofs, fixVals);
    % 
    %     Kuu_col_num = (RsolidPert - Rsolid0) / epsU;
    %     Kuu_col_ana = Kuu(:,j);
    % 
    %     err = norm(Kuu_col_ana - Kuu_col_num, inf) / ...
    %           max(norm(Kuu_col_num, inf), 1);
    % 
    %     fprintf('Kuu col %d rel error = %.3e, ||num||=%.3e, ||ana||=%.3e\n', ...
    %         j, err, norm(Kuu_col_num,inf), norm(Kuu_col_ana,inf));
    % end
    % % ============================================================
end

function merit = coupled_merit(solidNorm, fluidNorm, solidTarget, fluidTarget)
    merit = max(solidNorm / solidTarget, fluidNorm / fluidTarget);
end

function stop = stop_fsolve_when_coupled_converged( ...
    y, ~, state, old, meshE, interfaceE, z, par, ...
    free, fixDofs, fixVals, Ju, Jp, solidTarget, fluidTarget)

    stop = false;
    if ~strcmp(state, 'iter')
        return;
    end

    try
        [solidNorm, fluidNorm, gapMin] = coupled_residual_norms_from_y( ...
            y, old, meshE, interfaceE, z, par, free, fixDofs, fixVals, Ju, Jp);
        stop = solidNorm < solidTarget && fluidNorm < fluidTarget && gapMin > par.minGap;
        printFsolve = ~isfield(par, 'fsolveDisplay') || ~strcmpi(par.fsolveDisplay, 'off') || ...
            (isfield(par, 'debugVerbose') && par.debugVerbose);
        if stop && printFsolve
            fprintf(['   fsolve OutputFcn stop: solid=%.3e, fluid=%.3e, ', ...
                     'gap=%.3e\n'], solidNorm, fluidNorm, gapMin);
        end
    catch ME
        if ~(contains(ME.message, 'Negative or zero J') || ...
             contains(ME.message, 'Non-positive radius') || ...
             contains(ME.message, 'Element inverted') || ...
             contains(ME.message, 'Solid geometry guard failed') || ...
             contains(ME.message, 'Gap violates minGap'))
            rethrow(ME);
        end
    end
end

function check_coupled_jacobian(u, p, old, meshE, interfaceE, z, par, ...
    free, fixDofs, fixVals, Kuu, Kup, Kpu, Kpp, Rsolid0, Rfluid0)

    epsU = 1e-10;
    epsP = 1e-4;

    interfaceDofs = sort([2*interfaceE(:)-1; 2*interfaceE(:)]);
    [~, loc] = ismember(interfaceDofs, free);
    uCols = loc(loc > 0);
    if isempty(uCols)
        uCols = round(linspace(1, numel(free), min(8, numel(free))));
    else
        uCols = uCols(round(linspace(1, numel(uCols), min(8, numel(uCols)))));
    end

    fprintf('--- Coupled Jacobian check: displacement columns ---\n');
    for jj = 1:numel(uCols)
        col = uCols(jj);
        uPert = u;
        uPert(free(col)) = uPert(free(col)) + epsU;

        [RsPert, RfPert] = monolithic_analytical_residual_only( ...
            uPert, p, old, meshE, interfaceE, z, par, ...
            free, fixDofs, fixVals);

        numSolid = (RsPert - Rsolid0) / epsU;
        numFluid = (RfPert - Rfluid0) / epsU;

        relSolid = norm(Kuu(:,col) - numSolid, inf) / ...
            max([norm(numSolid, inf), norm(Kuu(:,col), inf), 1e-30]);
        relFluid = norm(Kpu(:,col) - numFluid, inf) / ...
            max([norm(numFluid, inf), norm(Kpu(:,col), inf), 1e-30]);

        fprintf('   u col %d: Kuu rel=%.3e, Kpu rel=%.3e\n', ...
                col, relSolid, relFluid);
    end

    pCols = unique(round(linspace(1, numel(p)-2, min(8, numel(p)-2))));

    fprintf('--- Coupled Jacobian check: pressure columns ---\n');
    for jj = 1:numel(pCols)
        col = pCols(jj);
        pPert = p;
        pPert(col+1) = pPert(col+1) + epsP;

        [RsPert, RfPert] = monolithic_analytical_residual_only( ...
            u, pPert, old, meshE, interfaceE, z, par, ...
            free, fixDofs, fixVals);

        numSolid = (RsPert - Rsolid0) / epsP;
        numFluid = (RfPert - Rfluid0) / epsP;

        relSolid = norm(Kup(:,col) - numSolid, inf) / ...
            max([norm(numSolid, inf), norm(Kup(:,col), inf), 1e-30]);
        relFluid = norm(Kpp(:,col) - numFluid, inf) / ...
            max([norm(numFluid, inf), norm(Kpp(:,col), inf), 1e-30]);

        fprintf('   p col %d: Kup rel=%.3e, Kpp rel=%.3e\n', ...
                col, relSolid, relFluid);
    end
end

function stateNew = build_coupled_state_from_up(u, p, old, meshE, interfaceE, z, par)
    N = numel(z);

    if isfield(par, 'limitSolidStep') && par.limitSolidStep
        duEq = u - old.uE;
        duMax = max(abs(duEq));
        alphaSolid = 1.0;

        if isfield(par, 'solidStepRelax') && isfinite(par.solidStepRelax)
            alphaSolid = min(alphaSolid, max(0.0, min(1.0, par.solidStepRelax)));
        end

        if isfield(par, 'maxSolidStepPerStep') && ...
                isfinite(par.maxSolidStepPerStep) && par.maxSolidStepPerStep > 0 && ...
                duMax > par.maxSolidStepPerStep
            alphaSolid = min(alphaSolid, par.maxSolidStepPerStep / duMax);
        end

        if alphaSolid < 1.0
            u = old.uE + alphaSolid * duEq;
            fprintf('   accepted solid relaxation alpha=%.3e, max|du|=%.3e m\n', ...
                alphaSolid, max(abs(u - old.uE)));
        end
    end

    stateNew = old;
    stateNew.uE = u;

    % Use the same endothelium interface geometry as the monolithic residual.
    [stateNew.deltaE, stateNew.UwE] = monolithic_interface_kinematics( ...
        meshE, u, old.uE, interfaceE, z, par);

    stateNew.deltaL = old.deltaL;
    stateNew.p = p;
    stateNew.pReduced = p;
    stateNew.uEPrev = old.uE;
    stateNew.pPrev = old.p;
    stateNew.dtPrev = par.dt;
    stateNew.UwL = par.UwL * ones(N,1);

    [~, ~, ~, tauL, tauE, uzL, uzE] = ...
        local_flux_and_shear(z, p, stateNew.deltaL, stateNew.deltaE, ...
                             stateNew.UwL, stateNew.UwE, par);

    stateNew.tauE = tauE;
    stateNew.tauL = tauL;
    stateNew.uzE = uzE;
    stateNew.uzL = uzL;

end

function [solidNorm, fluidNorm, gapMin] = coupled_residual_norms_from_y( ...
    y, old, meshE, interfaceE, z, par, free, fixDofs, fixVals, Ju, Jp)

    u = old.uE;
    p = old.p;

    u(free) = y(1:numel(free)) * Ju;
    u(fixDofs) = fixVals;

    p(2:end-1) = y(numel(free)+1:end) * Jp;
    p(1) = par.pIn;
    p(end) = par.pOut;

    [deltaE, ~] = monolithic_interface_kinematics_value_only( ...
        meshE, u, old.uE, interfaceE, z, par);
    gapMin = min(deltaE - old.deltaL);

    try
        assert_solid_geometry_ok( ...
            solid_geometry_quality(meshE, u, 'endothelium fsolve candidate'), par);
    catch
        solidNorm = inf;
        fluidNorm = inf;
        return;
    end

    if any(deltaE - old.deltaL <= par.minGap)
        solidNorm = inf;
        fluidNorm = inf;
        return;
    end

    [Rsolid, Rfluid] = monolithic_analytical_residual_only( ...
        u, p, old, meshE, interfaceE, z, par, free, fixDofs, fixVals);

    solidNorm = norm(Rsolid, inf);
    fluidNorm = norm(Rfluid, inf);
end

