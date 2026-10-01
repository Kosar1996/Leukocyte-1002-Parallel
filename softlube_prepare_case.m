function [par, S, uE_pre] = softlube_prepare_case(cfg)
%SOFTLUBE_PREPARE_CASE Load input geometry and build the solver parameter struct.
    if ~isfield(cfg, 'geometry') || ~isfield(cfg.geometry, 'endotheliumPrestressFile')
        error('cfg.geometry.endotheliumPrestressFile is required.');
    end

    S = load(cfg.geometry.endotheliumPrestressFile);
    par = S.par;
    par.endotheliumPrestressFile = cfg.geometry.endotheliumPrestressFile;
    uE_pre = S.uE_pre;

    par = softlube_original_default_parameters(par);
    par = softlube_apply_case_fields(par, cfg);
    par = softlube_finish_derived_parameters(par);

    useEnvironment = true;
    if isfield(cfg, 'runtime') && isfield(cfg.runtime, 'useEnvironment')
        useEnvironment = cfg.runtime.useEnvironment;
    end
    if useEnvironment
        par = softlube_apply_environment_overrides(par);
        par = softlube_finish_derived_parameters(par);
    end

    if isfield(cfg, 'parOverrides') && ~isempty(cfg.parOverrides)
        par = softlube_merge_struct(par, cfg.parOverrides);
        par = softlube_finish_derived_parameters(par);
    end

    fprintf(['FULL-DOF ANALYTICAL-JACOBIAN RUN: dt=%.3e, ', ...
             'full DOFs, adaptive retry %s.\n'], ...
             par.dt, string_on_off(par.enableAdaptiveTimeStep));
end

function par = softlube_original_default_parameters(par)

% Defaults from a2D_deformable_leukocyte_exact_interface.m before case overrides.
    par.usePrestressedLeukocyteIC = true;
    par.rigidLeukocyte = false;

    par.zMin = 0e-6;
    par.zMax = 4e-6;

    par.UwL = 0;
    par.UwE = 0;
    par.noLeukocyte = false;

    par.leukocytePrestressFile = 'solid_leu_P600.mat';
    par.useLeukocyteCenterline = true;
    par.RLin = 0;
    par.RLout = 4e-6;
    par.zMinL = -2e-6;
    par.zMaxL = 6e-6;

    par.dt = 1e-4;
    par.tEnd = 0.001;

    par.Ee = 500;
    par.nuE = 0.46;
    par.useViscoelasticEndothelium = true;
    par.etaE = 1;

    par.EL = 200;
    par.nuL = 0.46;
    par.useViscoelasticLeukocyte = true;
    par.etaL = 1;
    par.supportL = 'axis';

    par.useExactDeformedInterface = true;
    par.useExactDeformedInterfaceInMonolithic = true;
    par.useReferenceZForTractionMapping = false;
    par.zeroLeukocyteTractionOutsideOverlap = true;
    par.useObjectiveKelvinVoigt = true;
    par.useSlopeAwareWallKinematics = true;
    par.warnPrestressLoadMismatch = false;
    par.useFsolveOutputFcn = false;
    par.useDirectFluidSolve = true;

    par.mu = 1.2e-3;
    par.slipE = 0e-9;
    par.slipL = 0e-9;

    par.useFull2DFluid = true;
    par.NrFluid2D = 30;
    par.Nr = par.NrFluid2D;
    par.fluid2DPenaltyFactor = 1e6;
    par.fluid2DClampEnds = false;
    par.useHybridGap1DExterior2DFluid = false;
    par.hybridGapZ = [-0.2e-6, 4.2e-6];
    par.NrExterior2D = par.NrFluid2D;
    par.bodyFittedRadialFineWindow = [];
    par.bodyFittedRadialFineWeight = 4.0;
    par.hybridExteriorMinCells = 4;
    par.hybridExteriorFailMode = 'warn';
    par.hybridUseExterior2DPressureVector = true;
    par.useExactEndotheliumFluidDomain = false;

    par.useBodyFittedMACFluid = true;
    par.storeFull2DFluidHist = false;
    par.store2DFluidArrays = true;
    par.store2DPhysicalGridHist = true;
    par.store2DSpeedHist = true;
    par.useBodyFittedMACTractionCorrection = true;
    par.maxBodyFittedTractionCorrections = 1;
    par.bodyFittedTractionCorrectionRelax = 0.05;
    par.bodyFittedTractionCorrectionFailMode = 'warn';
    par.resolveFluidAfterTractionCorrection = true;
    par.useBodyFittedMACTractionInSolid = true;

    par.useRLoutFluidInterfaceForSolidLeukocyte = false;
    par.useRLoutFluidInterfaceForLeukocyte = false;
    par.useFixedCylindricalLeukocyte = false;

    % Default changed Sep 3 (Issue 2/3): solve_fluid_2D_bodyfitted_MAC.m
    % previously hardcoded its internal pressurePenalty to 0, silencing an
    % artificial-compressibility stabilization term already supported by
    % solve_stokes_bodyfitted_MAC.m. Directly measured (Sep 3) to improve
    % the fluid matrix's condition number 68x at the exact near-floor gap
    % state that used to cause the correction loop's persistent-oscillation
    % failure (85 passes, 18+ hours, never resolving). Validated across
    % three full scenarios (buggy and fixed rOuter, from normal early-run
    % conditions through the true 1nm gap floor): consistently prevents
    % that permanent-stall failure mode with no adverse effect observed in
    % any tested run. Still overridable via par.fluidPressurePenalty for
    % any config that needs a different value.
    par.fluidPressurePenalty = 10;

    par.minGap = 0.001e-6;
    par.relax0 = 0.001;
    par.relaxMin = 0.00001;
    par.relaxMax = 0.05;
    par.NrPlot = 1000;
    par.maxHist = 200;

    par.maxCoupling = 100;
    par.tolCoupling = 0.5e-4;
    par.maxLeukocyteFluidCoupling = 80;
    par.tolLeukocyteFluidCoupling = 5e-4;
    par.maxAcceptedLeukocyteCoupling = 1e-2;
    par.leukocyteCouplingRelax = 0.25;
    par.leukocyteCouplingRelaxMin = 0.02;
    par.maxNewtonFluid = 30;
    par.tolNewtonFluid = 1e-13;
    par.lineSearchMax = 80;
    par.SsrcFactor = 1;

    par.pIn = 0;
    par.pOut = 0;

    par.maxNewtonMono = 100;
    par.tolNewtonMono = 1e-6;
    par.lineSearchMaxMono = 80;
    par.fluidScaleFactorMono = 20;
    par.uScaleMono = 1e-6;
    par.pScaleMono = 100;
    par.trustU0 = 5e-8;
    par.trustUMin = 1e-12;
    par.trustUMax = 1e-6;

    if ~isfield(par, 'newtonMaxItSolid') || isempty(par.newtonMaxItSolid)
        par.newtonMaxItSolid = 300;
    else
        par.newtonMaxItSolid = max(par.newtonMaxItSolid, 300);
    end
    if ~isfield(par, 'newtonTolSolid') || isempty(par.newtonTolSolid)
        par.newtonTolSolid = 1e-6;
    end

    par.solidAbsTol = 1e-10;
    par.monoSolidAbsTol = 1e-7;
    par.monoFluidAbsTol = 1e-8;
    par.solidFallbackAbsTol = 1e-8;
    par.solidFallbackMinIterations = 2;
    par.solidTrustU0 = 2e-8;
    par.solidTrustUMin = 1e-13;
    par.solidTrustUMax = 2e-7;
    par.stopAtMinGap = true;
    par.gapStopFactor = 2;
    par.enableAdaptiveTimeStep = true;
    par.dtMin = par.dt / 1024;
    par.dtRetryFactor = 0.5;
    par.dtGrowFactor = 1.25;
    par.maxTimeStepRetries = 18;

    par.enablePressureJumpRetry = false;
    par.maxPressureJumpAbs = inf;
    par.maxPressureJump2DAbs = inf;
    par.maxPressureJumpRel = inf;
    par.pressureJumpUseMidpointOnly = false;
    par.pressureJumpMidHalfWidth = 2;
    par.pressureJumpIgnoreFirstStep = true;
    par.pressureJumpSafetyFactor = 1.05;
    par.acceptPressureJumpAtDtMin = true;

    par.minSolidJacobian = 1e-4;
    par.minSolidRadius = 1e-12;
    par.debugVerbose = false;
    par.debugCoupledJacobian = false;
    par.useFsolveMono = true;
    par.useMonoPredictor = false;
    par.fsolveAlgorithms = {'trust-region-dogleg', 'levenberg-marquardt'};
    par.fsolveDisplay = 'off';
    par.printEvery = 1;
    par.saveOutput = true;
    par.makePlots = true;
    par.plotFinalGlobalPressureContour = false;
    par.plotGlobal2DPressureTractionComparison = false;
    par.useGlobal2DPressureTraction = false;
    par.global2DPressureTractionBlend = 1.0;
    par.global2DPressureTractionNr = 81;
    par.global2DPressureTractionScreenLength = 0.5e-6;
    par.outputFile = 'simulation_output_two_solid_full_analytical_nopre.mat';
    % Periodic checkpoint save (Sep 11): independent of saveOutput/outputFile
    % above (which only write once, at the very end, after plotting) -- this
    % writes a partial 'out' during the run itself, so a wall-time kill or a
    % post-solve crash (e.g. in a plot call) doesn't lose an otherwise-good
    % run. Disabled by default (empty checkpointFile); case scripts opt in.
    par.checkpointFile = '';
    par.checkpointEvery = 10;

    par.twoSolidInterfaceOnlyTest = false;
    par.useTwoSolidSemiAnalyticalJacobian = true;
    par.useTwoSolidFullAnalyticalJacobian = true;
    par.checkTwoSolidAnalyticalJacobian = false;
    par.maxFunctionEvaluationsTwoSolid = 400;

    par.pressureLimiterEnabled = false;
    par.pressureLimiterMaxStepAbs = 25;
    par.diagnosticsEnabled = true;
    par.diagnosticsWarnFluidResidual = 1e-7;
    par.diagnosticsWarnVolumeJumpRel = inf;
    par.limitSolidStep = false;
    par.solidStepRelax = 0.5;
    par.maxSolidStepPerStep = 5e-8;
    par.NrL = 18;
end

function par = softlube_apply_case_fields(par, cfg)
    if isfield(cfg, 'geometry')
        g = cfg.geometry;
        par = softlube_copy_fields(par, g, { ...
            'leukocytePrestressFile', 'usePrestressedLeukocyteIC', ...
            'useLeukocyteCenterline', 'RLin', 'RLout', ...
            'zMinL', 'zMaxL', 'NzSolidL', 'NrL', 'Rc'});
    end

    if isfield(cfg, 'fluid')
        f = cfg.fluid;
        par = softlube_copy_fields(par, f, { ...
            'zMin', 'zMax', 'mu', 'slipE', 'slipL', ...
            'useFull2DFluid', 'NrFluid2D', 'Nr', ...
            'fluid2DPenaltyFactor', 'fluid2DClampEnds', ...
            'useHybridGap1DExterior2DFluid', 'hybridGapZ', ...
            'NrExterior2D', 'hybridExteriorMinCells', ...
            'bodyFittedRadialFineWindow', 'bodyFittedRadialFineWeight', ...
            'hybridExteriorFailMode', 'hybridUseExterior2DPressureVector', ...
            'useExactEndotheliumFluidDomain', ...
            'useBodyFittedMACFluid', 'storeFull2DFluidHist', ...
            'store2DFluidArrays', 'store2DPhysicalGridHist', ...
            'store2DSpeedHist', 'useBodyFittedMACTractionCorrection', ...
            'useFeedbackTractionCorrection', ...
            'debugRadialAxialTractionCorrection', ...
            'useAitkenTractionCorrectionRelax', ...
            'tractionCorrectionMismatchThresholdPct', ...
            'maxBodyFittedTractionCorrections', ...
            'bodyFittedTractionCorrectionRelax', ...
            'bodyFittedTractionCorrectionFailMode', ...
            'resolveFluidAfterTractionCorrection', ...
            'useBodyFittedMACTractionInSolid', ...
            'useRLoutFluidInterfaceForSolidLeukocyte', ...
            'useRLoutFluidInterfaceForLeukocyte', ...
            'useFixedCylindricalLeukocyte', ...
            'useSlopeAwareWallKinematics'});
        if isfield(f, 'innerBoundary')
            switch lower(f.innerBoundary)
                case {'deformed_leukocyte', 'state_deltal', 'exact'}
                    par.useFixedCylindricalLeukocyte = false;
                    par.useRLoutFluidInterfaceForSolidLeukocyte = false;
                    par.useRLoutFluidInterfaceForLeukocyte = false;
                case {'rlout', 'fixed_cylinder', 'reference_cylinder'}
                    par.useFixedCylindricalLeukocyte = true;
                    par.useRLoutFluidInterfaceForSolidLeukocyte = true;
                    par.useRLoutFluidInterfaceForLeukocyte = true;
                otherwise
                    error('Unknown cfg.fluid.innerBoundary: %s', f.innerBoundary);
            end
        end
    end

    if isfield(cfg, 'bc')
        par = softlube_copy_fields(par, cfg.bc, {'pIn', 'pOut', 'pInFun', 'pOutFun', 'UwL', 'UwE'});
    end

    if isfield(cfg, 'solid')
        s = cfg.solid;
        par = softlube_copy_fields(par, s, {'noLeukocyte', 'rigidLeukocyte'});
        if isfield(s, 'endothelium')
            par = softlube_copy_fields(par, s.endothelium, ...
                {'Ee', 'nuE', 'useViscoelasticEndothelium', 'etaE'});
        end
        if isfield(s, 'leukocyte')
            l = s.leukocyte;
            if isfield(l, 'enabled')
                par.noLeukocyte = ~l.enabled;
            end
            if isfield(l, 'deformable')
                par.rigidLeukocyte = ~l.deformable;
            end
            par = softlube_copy_fields(par, l, ...
                {'EL', 'nuL', 'useViscoelasticLeukocyte', 'etaL', 'supportL'});
        end
    end

    if isfield(cfg, 'interface')
        par = softlube_copy_fields(par, cfg.interface, { ...
            'useExactDeformedInterface', 'useExactDeformedInterfaceInMonolithic', ...
            'useReferenceZForTractionMapping', ...
            'zeroLeukocyteTractionOutsideOverlap', ...
            'leukocytePressureSupportZ', 'warnPrestressLoadMismatch'});
    end

    if isfield(cfg, 'numerics')
        par = softlube_copy_fields(par, cfg.numerics, { ...
            'dt', 'tEnd', 'minGap', 'relax0', 'relaxMin', 'relaxMax', ...
            'maxCoupling', 'tolCoupling', 'maxLeukocyteFluidCoupling', ...
            'tolLeukocyteFluidCoupling', 'maxAcceptedLeukocyteCoupling', ...
            'leukocyteCouplingRelax', 'leukocyteCouplingRelaxMin', ...
            'maxNewtonFluid', 'tolNewtonFluid', 'lineSearchMax', ...
            'SsrcFactor', 'maxNewtonMono', 'tolNewtonMono', ...
            'lineSearchMaxMono', 'fluidScaleFactorMono', ...
            'uScaleMono', 'pScaleMono', 'trustU0', 'trustUMin', ...
            'trustUMax', 'newtonMaxItSolid', 'newtonTolSolid', ...
            'solidAbsTol', 'monoSolidAbsTol', 'monoFluidAbsTol', ...
            'solidFallbackAbsTol', 'solidFallbackMinIterations', ...
            'solidTrustU0', 'solidTrustUMin', 'solidTrustUMax', ...
            'stopAtMinGap', 'gapStopFactor', 'enableAdaptiveTimeStep', ...
            'dtMin', 'dtRetryFactor', 'dtGrowFactor', ...
            'maxTimeStepRetries', 'enablePressureJumpRetry', ...
            'maxPressureJumpAbs', 'maxPressureJump2DAbs', ...
            'maxPressureJumpRel', 'pressureJumpUseMidpointOnly', ...
            'pressureJumpMidHalfWidth', 'pressureJumpIgnoreFirstStep', ...
            'pressureJumpSafetyFactor', 'acceptPressureJumpAtDtMin', ...
            'minSolidJacobian', 'minSolidRadius', 'debugVerbose', ...
            'debugCoupledJacobian', 'useFsolveMono', 'useMonoPredictor', ...
            'fsolveAlgorithms', 'fsolveDisplay', ...
            'twoSolidInterfaceOnlyTest', 'useTwoSolidSemiAnalyticalJacobian', ...
            'useTwoSolidFullAnalyticalJacobian', ...
            'checkTwoSolidAnalyticalJacobian', ...
            'maxFunctionEvaluationsTwoSolid', 'useObjectiveKelvinVoigt', ...
            'useFsolveOutputFcn', 'useDirectFluidSolve', ...
            'pressureLimiterEnabled', 'pressureLimiterMaxStepAbs', ...
            'diagnosticsEnabled', 'diagnosticsWarnFluidResidual', ...
            'diagnosticsWarnVolumeJumpRel', 'limitSolidStep', ...
            'solidStepRelax', 'maxSolidStepPerStep'});
    end

    if isfield(cfg, 'output')
        par = softlube_copy_fields(par, cfg.output, ...
            {'saveOutput', 'makePlots', 'printEvery', 'NrPlot', 'maxHist', ...
             'outputFile', 'storeFull2DFluidHist', 'store2DFluidArrays', ...
             'store2DPhysicalGridHist', 'store2DSpeedHist', ...
             'plotFinalGlobalPressureContour', ...
             'plotGlobal2DPressureTractionComparison', ...
             'checkpointFile', 'checkpointEvery'});
    end
end

function par = softlube_copy_fields(par, source, names)
    for k = 1:numel(names)
        name = names{k};
        if isfield(source, name)
            par.(name) = source.(name);
        end
    end
end

function par = softlube_finish_derived_parameters(par)
    par.Lz = par.zMax - par.zMin;
    if par.Lz <= 0
        error('Expected par.zMax > par.zMin.');
    end

    if isfield(par, 'zMinL') && isfield(par, 'zMaxL')
        par.LzL = par.zMaxL - par.zMinL;
    end
    if isfield(par, 'LzL') && isfield(par, 'NzSolid') && par.Lz > 0
        par.NzSolidL = round((par.NzSolid - 1) * par.LzL / par.Lz) + 1;
    end

    par.Ge = par.Ee/(2*(1+par.nuE));
    par.Ke = par.Ee/(3*(1-2*par.nuE));
    par.GL = par.EL/(2*(1+par.nuL));
    par.KL = par.EL/(3*(1-2*par.nuL));

    if isfield(par, 'useLeukocyteCenterline') && par.useLeukocyteCenterline
        par.RLin = 0;
        par.supportL = 'axis';
    end

    if isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure && ...
            isfield(par, 'global1DNzFluid') && isfinite(par.global1DNzFluid) && ...
            par.global1DNzFluid >= 3
        par.NzFluid = round(par.global1DNzFluid);
    else
        par.NzFluid = par.NzSolid;
    end
    par.dz = par.Lz/(par.NzFluid-1);

    if isfield(par, 'NrFluid2D') && (~isfield(par, 'Nr') || isempty(par.Nr))
        par.Nr = par.NrFluid2D;
    elseif isfield(par, 'Nr') && (~isfield(par, 'NrFluid2D') || isempty(par.NrFluid2D))
        par.NrFluid2D = par.Nr;
    end

    if isfield(par, 'NrFluid2D')
        par.Nr = par.NrFluid2D;
    end
    if isfield(par, 'fluid2DPenaltyFactor')
        par.penaltyLambda = par.fluid2DPenaltyFactor * par.mu / max(par.dt, realmin);
    end

    if use_hybrid_gap1d_exterior2d_fluid(par)
        par.useFull2DFluid = true;
        if ~isfield(par, 'NrExterior2D') || isempty(par.NrExterior2D)
            par.NrExterior2D = par.NrFluid2D;
        end
    end

    % ==================================== ADD THIS FIX ========================================
    % STRICT GUARD: Full 2D Fluid Mode explicitly overrides any residual hybrid flags
    if isfield(par, 'useFull2DFluid') && par.useFull2DFluid
        par.useHybridGap1DExterior2DFluid = false;
        par.useGlobal2DPressureTraction = false;
        par.useSmoothHybridBlending = false;
    end
% ==========================================================================================

    par.dtMin = min(par.dtMin, par.dt);
    if ~isfield(par, 'leukocytePressureSupportZ') || isempty(par.leukocytePressureSupportZ)
        par.leukocytePressureSupportZ = [par.zMin, par.zMax];
    end

    if par.useExactDeformedInterface && ...
            ~(isfield(par, 'useExactDeformedInterfaceInMonolithic') && ...
              par.useExactDeformedInterfaceInMonolithic)
        warning(['Exact deformed-interface post-update disabled because the ', ...
            'monolithic residual/Jacobian is linearized on the reference-z ', ...
            'interface. Enable only after adding exact-interface sensitivities.']);
        par.useExactDeformedInterface = false;
    end

    if isfield(par, 'useFull2DFluid') && par.useFull2DFluid && ...
            (~isfield(par, 'diagnosticsWarnFluidResidual') || ...
             par.diagnosticsWarnFluidResidual < 1e-7)
        par.diagnosticsWarnFluidResidual = 1e-7;
    end
end

function par = softlube_apply_environment_overrides(par)
    dtEnv = str2double(getenv('SOFTLUBE_DT'));
    if isfinite(dtEnv) && dtEnv > 0, par.dt = dtEnv; end

    tEndEnv = str2double(getenv('SOFTLUBE_TEND'));
    if isfinite(tEndEnv) && tEndEnv > 0, par.tEnd = tEndEnv; end

    full2DFluidEnv = strtrim(getenv('SOFTLUBE_FULL2D_FLUID'));
    if ~isempty(full2DFluidEnv)
        if any(strcmpi(full2DFluidEnv, {'1','true','yes','on','full2d'}))
            par.useFull2DFluid = true;
        elseif any(strcmpi(full2DFluidEnv, {'0','false','no','off','reynolds','lubrication'}))
            par.useFull2DFluid = false;
        end
    end

    nrFluid2DEnv = str2double(getenv('SOFTLUBE_NR_FLUID_2D'));
    if isfinite(nrFluid2DEnv) && nrFluid2DEnv >= 3
        par.NrFluid2D = round(nrFluid2DEnv);
        par.Nr = par.NrFluid2D;
    end

    par = softlube_env_bool(par, 'SOFTLUBE_PRESTRESSED_LEUKOCYTE', ...
        'usePrestressedLeukocyteIC', {'prestressed'}, {'unprestressed','undeformed'});

    leukocytePrestressFileEnv = strtrim(getenv('SOFTLUBE_LEUKOCYTE_PRESTRESS_FILE'));
    if ~isempty(leukocytePrestressFileEnv)
        par.leukocytePrestressFile = leukocytePrestressFileEnv;
    end

    par = softlube_env_bool(par, 'SOFTLUBE_RIGID_LEUKOCYTE', ...
        'rigidLeukocyte', {'rigid','fixed'}, {'deformable'});
    par = softlube_env_bool(par, 'SOFTLUBE_LEUKOCYTE_CENTERLINE', ...
        'useLeukocyteCenterline', {'axis','centerline'}, {'inner-radius'});
    par = softlube_env_bool(par, 'SOFTLUBE_ADAPTIVE_DT', ...
        'enableAdaptiveTimeStep', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_ENABLE_ADAPTIVE_DT', ...
        'enableAdaptiveTimeStep', {}, {});

    dtMinEnv = str2double(getenv('SOFTLUBE_DT_MIN'));
    if isfinite(dtMinEnv) && dtMinEnv > 0, par.dtMin = min(dtMinEnv, par.dt); end

    dtRetryFactorEnv = str2double(getenv('SOFTLUBE_DT_RETRY_FACTOR'));
    if isfinite(dtRetryFactorEnv) && dtRetryFactorEnv > 0 && dtRetryFactorEnv < 1
        par.dtRetryFactor = dtRetryFactorEnv;
    end

    dtGrowFactorEnv = str2double(getenv('SOFTLUBE_DT_GROW_FACTOR'));
    if isfinite(dtGrowFactorEnv) && dtGrowFactorEnv >= 1
        par.dtGrowFactor = dtGrowFactorEnv;
    end

    maxTimeStepRetriesEnv = str2double(getenv('SOFTLUBE_MAX_TIME_STEP_RETRIES'));
    if isfinite(maxTimeStepRetriesEnv) && maxTimeStepRetriesEnv >= 0
        par.maxTimeStepRetries = round(maxTimeStepRetriesEnv);
    end

    minSolidJacobianEnv = str2double(getenv('SOFTLUBE_MIN_SOLID_JACOBIAN'));
    if isfinite(minSolidJacobianEnv) && minSolidJacobianEnv >= 0
        par.minSolidJacobian = minSolidJacobianEnv;
    end

    minSolidRadiusEnv = str2double(getenv('SOFTLUBE_MIN_SOLID_RADIUS'));
    if isfinite(minSolidRadiusEnv) && minSolidRadiusEnv >= 0
        par.minSolidRadius = minSolidRadiusEnv;
    end

    monoSolidAbsTolEnv = str2double(getenv('SOFTLUBE_MONO_SOLID_ABS_TOL'));
    if isfinite(monoSolidAbsTolEnv) && monoSolidAbsTolEnv > 0
        par.monoSolidAbsTol = monoSolidAbsTolEnv;
    end

    monoFluidAbsTolEnv = str2double(getenv('SOFTLUBE_MONO_FLUID_ABS_TOL'));
    if isfinite(monoFluidAbsTolEnv) && monoFluidAbsTolEnv > 0
        par.monoFluidAbsTol = monoFluidAbsTolEnv;
    end

    printEveryEnv = str2double(getenv('SOFTLUBE_PRINT_EVERY'));
    if isfinite(printEveryEnv) && printEveryEnv > 0
        par.printEvery = max(1, round(printEveryEnv));
    end

    fsolveDisplayEnv = strtrim(getenv('SOFTLUBE_FSOLVE_DISPLAY'));
    if ~isempty(fsolveDisplayEnv), par.fsolveDisplay = fsolveDisplayEnv; end

    par = softlube_env_bool(par, 'SOFTLUBE_SAVE_OUTPUT', 'saveOutput', {}, {});
    skipSaveEnv = strtrim(getenv('SOFTLUBE_SKIP_SAVE'));
    if any(strcmpi(skipSaveEnv, {'1','true','yes','on'})), par.saveOutput = false; end

    par = softlube_env_bool(par, 'SOFTLUBE_MAKE_PLOTS', 'makePlots', {}, {});
    skipPlotsEnv = strtrim(getenv('SOFTLUBE_SKIP_PLOTS'));
    if any(strcmpi(skipPlotsEnv, {'1','true','yes','on'})), par.makePlots = false; end

    par = softlube_env_bool(par, 'SOFTLUBE_PRESSURE_LIMITER', ...
        'pressureLimiterEnabled', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_DIAGNOSTICS', ...
        'diagnosticsEnabled', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_EXACT_DEFORMED_INTERFACE', ...
        'useExactDeformedInterface', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_EXACT_ENDOTHELIUM_DOMAIN', ...
        'useExactEndotheliumFluidDomain', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_OBJECTIVE_KV', ...
        'useObjectiveKelvinVoigt', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_DIRECT_FLUID', ...
        'useDirectFluidSolve', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_USE_FSOLVE_MONO', ...
        'useFsolveMono', {}, {});
    par = softlube_env_bool(par, 'SOFTLUBE_PRESSURE_JUMP_RETRY', ...
        'enablePressureJumpRetry', {}, {});

    debugVerboseEnv = strtrim(getenv('SOFTLUBE_DEBUG_VERBOSE'));
    if any(strcmpi(debugVerboseEnv, {'1','true','yes','on'}))
        par.debugVerbose = true;
    end

    fastModeEnv = strtrim(getenv('SOFTLUBE_FAST'));
    if any(strcmpi(fastModeEnv, {'1','true','yes','on'}))
        par.NrL = min(par.NrL, 30);
        par.maxNewtonMono = min(par.maxNewtonMono, 80);
        par.tolNewtonFluid = max(par.tolNewtonFluid, 1e-11);
        par.solidAbsTol = max(par.solidAbsTol, 1e-11);
        par.solidFallbackAbsTol = max(par.solidFallbackAbsTol, 1e-8);
        par.newtonTolSolid = max(par.newtonTolSolid, 1e-6);
        par.newtonMaxItSolid = min(par.newtonMaxItSolid, 40);
        par.maxLeukocyteFluidCoupling = min(par.maxLeukocyteFluidCoupling, 8);
        par.tolLeukocyteFluidCoupling = max(par.tolLeukocyteFluidCoupling, 2e-3);
        par.maxAcceptedLeukocyteCoupling = max(par.maxAcceptedLeukocyteCoupling, 2e-2);
    end
end

function par = softlube_env_bool(par, envName, fieldName, trueAliases, falseAliases)
    value = strtrim(getenv(envName));
    if isempty(value), return; end
    trueWords = [{'1','true','yes','on'}, trueAliases(:).'];
    falseWords = [{'0','false','no','off'}, falseAliases(:).'];
    if any(strcmpi(value, trueWords))
        par.(fieldName) = true;
    elseif any(strcmpi(value, falseWords))
        par.(fieldName) = false;
    end
end

function out = softlube_merge_struct(out, overrides)
    names = fieldnames(overrides);
    for k = 1:numel(names)
        out.(names{k}) = overrides.(names{k});
    end
end

