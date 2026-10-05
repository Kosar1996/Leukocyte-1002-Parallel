%% =========================================================================
% HEADER SUMMARY OF CHANGES:
% 1. Section 4 (Fluid Mesh Resolution): Expanded the gap1DWindow boundaries 
%    from [-0.2e-6, 4.2e-6] to [-0.5e-6, 4.5e-6] to create a smoother spatial 
%    buffer and eliminate sharp velocity profile discontinuities at t = 0.
% 2. Section 8 (Hybrid Gap-1D / Exterior-2D Fluid Module): Introduced smooth 
%    hybrid blending parameters (`hybridTransitionBuffer` and `useSmoothHybridBlending`) 
%    to grade the transition between the 1D lubrication and full 2D solvers.
% =========================================================================

%% RUN_TEMPLATE
% Base case per item 1's validated codebase + item 1's original meshfiles
% (solid_leukocyte_P600.mat, solid_endothelium_P300.mat), 1000 steps or
% 4 hours, whichever first.[cite: 5]

clear functions;
clear classes;
clear all;
rehash toolbox;
rehash path;
java.lang.System.gc();
addpath(pwd, '-begin');
clc;

fprintf('\n=== MATLAB Function Lookup Diagnostics ===\n');
which('-all', 'solve_finite_def_solid');
which('-all', 'apply_bodyfitted_MAC_traction_correction');
fprintf('===========================================\n\n');




%% 1. Paths
softlubeDir = '';
% FIX 1: Resolve softlubeDir dynamically relative to mfilename if empty
if isempty(softlubeDir)
    softlubeDir = fileparts(mfilename('fullpath'));
end
if ~isempty(softlubeDir)
    addpath(softlubeDir);
end

% INPUT 1: Path to the MAT file containing historical step data
matFilePath = fullfile(softlubeDir, 'case_7_t_step14.mat');

% INPUT 2: Exact timestep index to roll back to and resume from
targetStep = 13;

%% 2. CUSTOMIZABLE PARAMETERS
dt = 2.5e-6;
nSteps = 1000;

leukocyteDeformable = true;
leukocyteUseViscoelastic = true;
 
if ~isempty(softlubeDir)
    addpath(softlubeDir);
end

% INPUT 1: Path to the MAT file containing historical step data
matFilePath = fullfile(softlubeDir, 'case_7_t_step14.mat');

% INPUT 2: Exact timestep index to roll back to and resume from
targetStep = 13;

%% 2. CUSTOMIZABLE PARAMETERS
dt = 2.5e-6;
nSteps = 1000;

leukocyteDeformable = true;
leukocyteUseViscoelastic = true;
leukocyteEtaL = 1;

% Late-Stage TEM Active Translocation Parameters
useActiveTranslocation = true;
fz_active_translocation = -1.0e9; % Active force density [N/m^3]
z_pore_bottom = 3.5e-6;          % z-coordinate of bottom exit [m]
sigma_head_decay = 0.5e-6;  
leukocyteEtaL = 1;


endotheliumUseViscoelastic = true;
endotheliumEtaE = 1;

tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 100000;

%% 3. Global Domain and Pressure Boundary Conditions
zMin = -6e-6;
zMax = 10e-6;
rOuter = 6e-6;

pAtZMin = 0;
pAtZMax = 0;
axisPressureBC = 'dPdr=0';

%% 4. Fluid Mesh Resolution
% OLD: gap1DWindow = [-0.2e-6, 4.2e-6];[cite: 5]
% FIX: Expanded window to provide a smooth buffer region preventing sharp velocity jumps
gap1DWindow = [-0.5e-6, 4.5e-6];
finePressureWindow = gap1DWindow;
coarseDz = 1.0e-6;
fineDz = 0.05e-6;

NrExterior2D = 20;
NrFluid2D = NrExterior2D;
radialFineWindow = [2e-6, rOuter];
radialFineWeight = 4.0;

%% 5. Plotting and Output
plotNative2DPressure = true;
plotNative2DStress = true;
plotNative2DVelocity = false;
plotFirstStepHybridMesh = false;
plotGlobalDomainSchematic = false;
makeLegacyPlots = false;
saveOutput = true;
outputFile = fullfile(softlubeDir, 'case_7_t.mat');
closeFiguresAtStart = true;
printEvery = 1;

%% 6. Build the Base Case
cfg = struct();

cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endothelium_P300.mat');
cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leukocyte_P300.mat');
cfg.geometry.usePrestressedLeukocyteIC = true;
cfg.geometry.useLeukocyteCenterline = true;

cfg.fluid.zMin = zMin;
cfg.fluid.zMax = zMax;
cfg.fluid.mu = 1.2e-3;
cfg.fluid.slipE = 0e-9;
cfg.fluid.slipL = 0e-9;
cfg.fluid.NrFluid2D = NrFluid2D;
cfg.fluid.Nr = NrFluid2D;
cfg.fluid.NrExterior2D = NrExterior2D;
cfg.fluid.bodyFittedRadialFineWindow = radialFineWindow;
cfg.fluid.bodyFittedRadialFineWeight = radialFineWeight;
cfg.fluid.innerBoundary = 'deformed_leukocyte';
cfg.fluid.useSlopeAwareWallKinematics = true;
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.2; % Reduced from 0.5 for smoother ICs

cfg.bc.pIn = pAtZMin;
cfg.bc.pOut = pAtZMax;
cfg.bc.UwL = 0;
cfg.bc.UwE = 0;

cfg.solid.noLeukocyte = false;
cfg.solid.leukocyte.enabled = true;
cfg.solid.leukocyte.deformable = leukocyteDeformable;
cfg.solid.leukocyte.useViscoelasticLeukocyte = leukocyteUseViscoelastic;
cfg.solid.leukocyte.etaL = leukocyteEtaL;
cfg.solid.leukocyte.supportL = 'axis';

% Late-Stage TEM Active Translocation Parameters
cfg.solid.leukocyte.useActiveTranslocation  = useActiveTranslocation;
cfg.solid.leukocyte.fz_active_translocation = fz_active_translocation;
cfg.solid.leukocyte.z_pore_bottom           = z_pore_bottom;
cfg.solid.leukocyte.sigma_head_decay        = sigma_head_decay;

cfg.solid.endothelium.useViscoelasticEndothelium = endotheliumUseViscoelastic;
cfg.solid.endothelium.etaE = endotheliumEtaE;

% cfg.interface.useExactDeformedInterface = true;
% cfg.interface.useExactDeformedInterfaceInMonolithic = true;
% cfg.interface.useReferenceZForTractionMapping = false;
% cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
% cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
% cfg.interface.warnPrestressLoadMismatch = false;
% cfg.numerics.dt = dt;
% cfg.numerics.tEnd = tEnd;
% cfg.numerics.minGap = 0.001e-6;
% cfg.numerics.maxNewtonMono = 100;
% cfg.numerics.tolNewtonMono = 1e-6;
% cfg.numerics.monoSolidAbsTol = 1e-7;
% cfg.numerics.monoFluidAbsTol = 1e-8;
% cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
% cfg.numerics.enablePressureJumpRetry = false;
% cfg.numerics.maxPressureJumpAbs = inf;
% cfg.numerics.maxPressureJump2DAbs = inf;
% cfg.numerics.useObjectiveKelvinVoigt = true;
% cfg.numerics.useDirectFluidSolve = true;
% cfg.numerics.useFsolveMono = true;
% cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
% cfg.numerics.dtMin = dtMin;
% cfg.numerics.dtRetryFactor = 0.5;
% cfg.numerics.dtGrowFactor = 1.25;
% cfg.numerics.maxTimeStepRetries = 18;
% cfg.numerics.fsolveDisplay = 'off';
% cfg.numerics.diagnosticsEnabled = true;

% ========================= NEW CODE ===========================
cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;
cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
cfg.interface.warnPrestressLoadMismatch = false;

% STEP 2 FIX: Smooth traction interpolation across interface boundaries
cfg.interface.tractionInterpMethod = 'pchip'; 

cfg.numerics.dt = dt;
cfg.numerics.tEnd = tEnd;
cfg.numerics.minGap = 0.001e-6;
cfg.numerics.maxNewtonMono = 300; % Increased from 100
cfg.numerics.tolNewtonMono = 1e-6;
cfg.numerics.monoSolidAbsTol = 1e-7;
cfg.numerics.monoFluidAbsTol = 1e-8;
cfg.numerics.maxFunctionEvaluationsTwoSolid = 600;
cfg.numerics.enablePressureJumpRetry = true;
cfg.numerics.maxPressureJumpAbs = 2.0e3;
cfg.numerics.maxPressureJump2DAbs = inf;
cfg.numerics.useObjectiveKelvinVoigt = true;
cfg.numerics.useDirectFluidSolve = true;
cfg.numerics.useFsolveMono = true;
cfg.numerics.enableICRamp = true;
cfg.numerics.icRampSteps = 10; % Smoothly ramps solid kinematic BCs over steps 1-10
cfg.numerics.limitSolidStep = true;
cfg.numerics.maxSolidStepPerStep = 1.0e-8; % 10 nm max nodal displacement per step
cfg.numerics.solidStepRelax = 0.5;

% STEP 1 FIX: Solver Damping & Adaptive Guard
cfg.numerics.fsolveDampingFactor = 0.5; % Damps step size to prevent overshooting
cfg.numerics.maxLineSearchStep = 0.25;   % Caps nodal displacement per sub-iteration
cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
cfg.numerics.dtMin = dtMin;
cfg.numerics.dtRetryFactor = 0.5;
cfg.numerics.dtGrowFactor = 1.1;        % Controlled step growth (reduced from 1.25)
cfg.numerics.maxTimeStepRetries = 18;
cfg.numerics.fsolveDisplay = 'iter';
cfg.numerics.diagnosticsEnabled = true;
% ==============================================================

cfg.output.saveOutput = saveOutput;
cfg.output.outputFile = outputFile;
cfg.output.checkpointFile = outputFile;
cfg.output.checkpointEvery = 1;
cfg.output.makePlots = makeLegacyPlots;
cfg.output.plotFinalGlobalPressureContour = false;
cfg.output.plotGlobal2DPressureTractionComparison = false;
cfg.output.printEvery = printEvery;
cfg.output.storeFull2DFluidHist = false;
cfg.output.store2DFluidArrays = true;
cfg.output.store2DPhysicalGridHist = true;
cfg.output.store2DSpeedHist = true;

cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = closeFiguresAtStart;
cfg.parOverrides = struct();
% ========================= NEW CODE ===========================
% STRATEGY A MASTER SAFEGUARDS & NUMERICAL TOLERANCES
cfg.numerics.maxFluidPressureCap = 3000.0; % [Pa] Upper bound on fluid pressure
cfg.numerics.maxFluidShearCap    = 500.0;  % [Pa] Upper bound on wall shear stress
cfg.numerics.v_max_cap           = 0.020;   % [m/s] Kinematic velocity clamp (20 mm/s)
cfg.numerics.solidAbsTol         = 1.0e-12; % [N] Absolute force tolerance
cfg.numerics.solidFallbackAbsTol = 5.0e-9;  % [N] Micro-element floor (5 nN)
cfg.numerics.newtonTolSolid      = 1.0e-3;  % Relative Newton force tolerance
cfg.numerics.eta_solid           = 10.0;    % [Pa*s] Viscoelastic rate dissipation

% STRATEGY A OUTER COUPLING & PREDICTOR CONTROLS
cfg.numerics.omega_min  = 0.08;   % Lower floor prevents lock-in
cfg.numerics.omega_max  = 0.50;   % Upper stability clamp
cfg.numerics.omega_init = 0.08;   % Initial Aitken factor
cfg.numerics.alpha_pred = 0.30;   % Damped warm-start predictor factor
% ==============================================================

% Sync Strategy A parameters directly into parOverrides
cfg.parOverrides.omega_min  = cfg.numerics.omega_min;
cfg.parOverrides.omega_max  = cfg.numerics.omega_max;
cfg.parOverrides.omega_init = cfg.numerics.omega_init;
cfg.parOverrides.alpha_pred = cfg.numerics.alpha_pred;

%% 7. Global Axial Mesh and Exterior Boundary Embedding
cfg.global1D.zDomain = [zMin, zMax];
cfg.global1D.rDomain = [0, rOuter];
cfg.global1D.axisPressureBC = axisPressureBC;
cfg.global1D.outerPressureBC = 0;

cfg.parOverrides.useGlobal1DPressure = true;
cfg.parOverrides.global1DOuterRadius = rOuter;
cfg.parOverrides.global1DAxisRadius = 1e-9;
cfg.parOverrides.useGlobal1DCoarseEdgeMesh = true;
cfg.parOverrides.global1DFineWindow = finePressureWindow;
cfg.parOverrides.global1DCoarseDz = coarseDz;
cfg.parOverrides.global1DFineDz = fineDz;
cfg.parOverrides.global1DPlotNr = 161;
cfg.parOverrides.fsolveDisplay = 'iter';

%% 8. Pure Full 2D Fluid Module Configuration
cfg.fluid.useFull2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = false;
cfg.fluid.useSmoothHybridBlending = false;

cfg.fluid.useBodyFittedMACFluid = true;
cfg.fluid.useBodyFittedMACTractionCorrection = true;
cfg.fluid.useBodyFittedMACTractionInSolid = true;
cfg.fluid.maxBodyFittedTractionCorrections = 20;
cfg.fluid.bodyFittedTractionCorrectionTol = 1e-2; % Increase tolerance from 1e-4 during gap squeeze
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.08; % 0.8% base relaxation for high-pressure squeezing
cfg.fluid.useAitkenTractionCorrectionRelax = true; 
cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.fluid.resolveFluidAfterTractionCorrection = true;
cfg.fluid.maxTractionNormChangeRatio = 2.0; % Rejects sub-pass if traction surges > 2x in 1 iteration

% Synchronize cfg.parOverrides directly with cfg.fluid & cfg.numerics
cfg.parOverrides.useFull2DFluid = true;
cfg.parOverrides.useHybridGap1DExterior2DFluid = false;
cfg.parOverrides.useGlobal2DPressureTraction = false; % CRITICAL: Prevent hybrid mode trigger

cfg.parOverrides.useBodyFittedMACFluid = true;
cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
cfg.parOverrides.maxBodyFittedTractionCorrections = 20;
cfg.parOverrides.bodyFittedTractionCorrectionTol = 1e-2; % Increase tolerance from 1e-4 during gap squeeze
cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.08;
cfg.parOverrides.useAitkenTractionCorrectionRelax = true;
cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
cfg.parOverrides.maxTractionNormChangeRatio = 2.0;

cfg.parOverrides.NrFluid2D = NrFluid2D;
cfg.parOverrides.Nr = NrFluid2D;
cfg.parOverrides.NrExterior2D = NrExterior2D;
cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
cfg.parOverrides.maxPressureJumpAbs = 2.0e3;
cfg.parOverrides.monoFluidAbsTol = 1e-8;

% Solver & Time-stepping Overrides (FIXED: Pass numerics settings into parOverrides)
cfg.parOverrides.tractionInterpMethod = cfg.interface.tractionInterpMethod;
cfg.parOverrides.fsolveDampingFactor = cfg.numerics.fsolveDampingFactor;
cfg.parOverrides.maxLineSearchStep = cfg.numerics.maxLineSearchStep;
cfg.parOverrides.maxNewtonMono = cfg.numerics.maxNewtonMono; % Passes 300 iteration cap
cfg.parOverrides.dtMin = dtMin;                             % Passes 1e-10 s minimum dt

% Solid Step Limiting
cfg.parOverrides.limitSolidStep = true;
cfg.parOverrides.maxSolidStepPerStep = 1.0e-8; % 10 nm cap
cfg.parOverrides.solidStepRelax = 0.2;

% Synchronize directly into parOverrides for restart mode
cfg.parOverrides.checkpointFile = outputFile;
cfg.parOverrides.checkpointEvery = 1;
cfg.parOverrides.saveOutput = saveOutput;
cfg.parOverrides.outputFile = outputFile;
% =========================================

% ========================= NEW CODE ===========================
% Synchronize ALL Strategy A Parameters into cfg.parOverrides
cfg.parOverrides.maxFluidPressureCap = cfg.numerics.maxFluidPressureCap;
cfg.parOverrides.maxFluidShearCap    = cfg.numerics.maxFluidShearCap;
cfg.parOverrides.v_max_cap           = cfg.numerics.v_max_cap;
cfg.parOverrides.solidAbsTol         = cfg.numerics.solidAbsTol;
cfg.parOverrides.solidFallbackAbsTol = cfg.numerics.solidFallbackAbsTol;
cfg.parOverrides.newtonTolSolid      = cfg.numerics.newtonTolSolid;
cfg.parOverrides.eta_solid           = cfg.numerics.eta_solid;
% In Section 6 & 8 of run_0928_2.m:
cfg.parOverrides.useActiveTranslocation  = useActiveTranslocation;
cfg.parOverrides.fz_active_translocation = fz_active_translocation;
cfg.parOverrides.z_pore_bottom           = z_pore_bottom;
cfg.parOverrides.sigma_head_decay        = sigma_head_decay;

% ==============================================================

%% 9. Runtime Execution & Rollback Launch
if plotNative2DPressure || plotNative2DVelocity || ...
        plotFirstStepHybridMesh || plotGlobalDomainSchematic
    set(0, 'DefaultFigureVisible', 'on');
end



fprintf('\nRunning 06_run_file_fixed: Strategy A Safeguarded Codebase\n');
fprintf('   dt                  = %.6e s\n', dt);
fprintf('   tEnd                = %.6e s\n', tEnd);
fprintf('   requested steps     = %.0f\n', nSteps);

% Print Strategy A verification check
fprintf('\n=== RUNTIME STRATEGY A CONFIGURATION CHECK ===\n');
fprintf('   omega_min  = %.2f (elevated floor)\n', cfg.numerics.omega_min);
fprintf('   omega_max  = %.2f (upper clamp)\n', cfg.numerics.omega_max);
fprintf('   alpha_pred = %.2f (damped predictor)\n', cfg.numerics.alpha_pred);
fprintf('   v_max_cap  = %.1e m/s (kinematic ceiling)\n', cfg.numerics.v_max_cap);
fprintf('===============================================\n\n');

% Add right before calling softlube_run_case_global_coupled:
fprintf('\n=========================================================\n');
fprintf('  RESTART EXECUTION DIAGNOSTICS\n');
fprintf('  - Target Restart Step : %d\n', targetStep);
fprintf('  - Input File Path     : %s\n', matFilePath);
fprintf('  - Traction Caps       : Pressure <= %.1f Pa | Shear <= %.1f Pa\n', ...
    cfg.numerics.maxFluidPressureCap, cfg.numerics.maxFluidShearCap);
fprintf('  - Pressure Cap        : maxP <= %.1f Pa\n', cfg.numerics.maxFluidPressureCap);
fprintf('  - Kinematic Clamp     : Velocity <= %.1f mm/s\n', cfg.numerics.v_max_cap * 1e3);
fprintf('  - Solid Viscosity     : eta_solid = %.1f Pa*s\n', cfg.numerics.eta_solid);
fprintf('  - Active Propulsion   : %s\n', char(string(cfg.solid.leukocyte.useActiveTranslocation)));
fprintf('  - Active Force f0     : %.3e N/m^3\n', cfg.solid.leukocyte.fz_active_translocation);
fprintf('  - Pore Exit z_bot     : %.3e m\n', cfg.solid.leukocyte.z_pore_bottom);
fprintf('  - Decay Length sigma  : %.3e m\n', cfg.solid.leukocyte.sigma_head_decay);
fprintf('=========================================================\n\n');

if exist(matFilePath, 'file')
    fprintf('[RESTART] Executing softlube_run_case_global_coupled(''%s'', %d, cfg)...\n', ...
        matFilePath, targetStep);
    % Pass runtime cfg struct as 3rd argument so overrides are applied on restart
    out = softlube_run_case_global_coupled(matFilePath, targetStep, cfg);
else
    fprintf('[FRESH START] Launching from t = 0 s using cfg...\n');
    out = softlube_run_case_global_coupled(cfg);
end

out.cfg = cfg;
save(outputFile, 'out', '-v7.3');

try
    out.native2D = collect_final_native2d(out);
    if plotNative2DPressure
        plot_select_native2d_pressure(out,out.stopStep);
    end
    if plotNative2DStress
        plot_select_native2d_stress(out,out.stopStep);
    end

    fprintf('\n06_run_file_fixed run finished\n');
    fprintf('   stopStep = %d\n', out.stopStep);
    fprintf('   final t  = %.6e s\n', out.t(end));
    if isfield(out, 'p2DMaxHist') && any(isfinite(out.p2DMaxHist))
        fprintf('   final max |P2D| = %.6e Pa\n', out.p2DMaxHist(end));
    end
catch MEpost
    warning(['Post-processing (plots/summary) failed, but the simulation ' ...
        'data was already saved above and is not affected: %s'], MEpost.message);
end

function native2D = collect_final_native2d(out)
native2D = struct();
if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
        ~isfield(out, 'PHist') || isempty(out.PHist)
    return;
end
k = out.stopStep;
native2D.P = out.PHist(:,:,k);
native2D.R = out.RPHist(:,:,k);
native2D.Z = out.ZPHist(:,:,k);
native2D.ur = out.urCHist(:,:,k);
native2D.uz = out.uzCHist(:,:,k);
native2D.speed = out.speedCHist(:,:,k);
end
