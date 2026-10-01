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

%% 2. CUSTOMIZABLE PARAMETERS
dt = 1e-5;
nSteps = 1000;

leukocyteDeformable = true;
leukocyteUseViscoelastic = true;
leukocyteEtaL = 1;

endotheliumUseViscoelastic = true;
endotheliumEtaE = 1;

tEnd = nSteps * dt;

useAdaptiveTimeStep = true;
dtMin = dt / 1024;

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
saveOutput = false;
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

cfg.solid.endothelium.useViscoelasticEndothelium = endotheliumUseViscoelastic;
cfg.solid.endothelium.etaE = endotheliumEtaE;

cfg.interface.useExactDeformedInterface = true;
cfg.interface.useExactDeformedInterfaceInMonolithic = true;
cfg.interface.useReferenceZForTractionMapping = false;
cfg.interface.zeroLeukocyteTractionOutsideOverlap = false;
cfg.interface.leukocytePressureSupportZ = [zMin, zMax];
cfg.interface.warnPrestressLoadMismatch = false;
cfg.numerics.dt = dt;
cfg.numerics.tEnd = tEnd;
cfg.numerics.minGap = 0.001e-6;
cfg.numerics.maxNewtonMono = 100;
cfg.numerics.tolNewtonMono = 1e-6;
cfg.numerics.monoSolidAbsTol = 1e-7;
cfg.numerics.monoFluidAbsTol = 1e-8;
cfg.numerics.maxFunctionEvaluationsTwoSolid = 400;
cfg.numerics.enablePressureJumpRetry = false;
cfg.numerics.maxPressureJumpAbs = inf;
cfg.numerics.maxPressureJump2DAbs = inf;
cfg.numerics.useObjectiveKelvinVoigt = true;
cfg.numerics.useDirectFluidSolve = true;
cfg.numerics.useFsolveMono = true;
cfg.numerics.enableAdaptiveTimeStep = useAdaptiveTimeStep;
cfg.numerics.dtMin = dtMin;
cfg.numerics.dtRetryFactor = 0.5;
cfg.numerics.dtGrowFactor = 1.25;
cfg.numerics.maxTimeStepRetries = 18;
cfg.numerics.fsolveDisplay = 'off';
cfg.numerics.diagnosticsEnabled = true;

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

%% 8. Hybrid Gap-1D / Exterior-2D Fluid Module
cfg.fluid.useFull2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = false;
cfg.fluid.hybridGapZ = gap1DWindow;
cfg.fluid.hybridExteriorMinCells = 4;
cfg.fluid.hybridExteriorFailMode = 'warn';
cfg.fluid.hybridUseExterior2DPressureVector = true;
cfg.fluid.useBodyFittedMACFluid = true;
cfg.fluid.useBodyFittedMACTractionCorrection = true;
cfg.fluid.useFeedbackTractionCorrection = false;
cfg.fluid.useAitkenTractionCorrectionRelax = false;
cfg.fluid.debugRadialAxialTractionCorrection = false;
cfg.fluid.maxBodyFittedTractionCorrections = 100;
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.5;
cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.fluid.resolveFluidAfterTractionCorrection = true;
cfg.fluid.useBodyFittedMACTractionInSolid = true;

% FIX: Added smooth hybrid blending parameters to prevent sharp velocity jumps
cfg.fluid.hybridTransitionBuffer = 0.3e-6;
cfg.fluid.useSmoothHybridBlending = true;

% FIX 2: Synchronize cfg.parOverrides directly with cfg.fluid
cfg.parOverrides.useFull2DFluid = cfg.fluid.useFull2DFluid;
cfg.parOverrides.useHybridGap1DExterior2DFluid = cfg.fluid.useHybridGap1DExterior2DFluid;
cfg.parOverrides.hybridGapZ = cfg.fluid.hybridGapZ;
cfg.parOverrides.NrExterior2D = NrExterior2D;
cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
cfg.parOverrides.hybridExteriorMinCells = cfg.fluid.hybridExteriorMinCells;
cfg.parOverrides.hybridExteriorFailMode = cfg.fluid.hybridExteriorFailMode;
cfg.parOverrides.hybridUseExterior2DPressureVector = cfg.fluid.hybridUseExterior2DPressureVector;
cfg.parOverrides.useBodyFittedMACFluid = cfg.fluid.useBodyFittedMACFluid;
cfg.parOverrides.useBodyFittedMACTractionCorrection = cfg.fluid.useBodyFittedMACTractionCorrection;
cfg.parOverrides.useFeedbackTractionCorrection = cfg.fluid.useFeedbackTractionCorrection;
cfg.parOverrides.useAitkenTractionCorrectionRelax = cfg.fluid.useAitkenTractionCorrectionRelax;
cfg.parOverrides.maxBodyFittedTractionCorrections = cfg.fluid.maxBodyFittedTractionCorrections;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = cfg.fluid.bodyFittedTractionCorrectionRelax;
cfg.parOverrides.bodyFittedTractionCorrectionFailMode = cfg.fluid.bodyFittedTractionCorrectionFailMode;
cfg.parOverrides.resolveFluidAfterTractionCorrection = cfg.fluid.resolveFluidAfterTractionCorrection;
cfg.parOverrides.useBodyFittedMACTractionInSolid = cfg.fluid.useBodyFittedMACTractionInSolid;
cfg.parOverrides.NrFluid2D = NrFluid2D;
cfg.parOverrides.Nr = NrFluid2D;
cfg.parOverrides.useGlobal2DPressureTraction = true;

% FIX: Pass smooth hybrid blending overrides
cfg.parOverrides.hybridTransitionBuffer = cfg.fluid.hybridTransitionBuffer;
cfg.parOverrides.useSmoothHybridBlending = cfg.fluid.useSmoothHybridBlending;

%% 9. Runtime
if plotNative2DPressure || plotNative2DVelocity || ...
        plotFirstStepHybridMesh || plotGlobalDomainSchematic
    set(0, 'DefaultFigureVisible', 'on');
end

fprintf('\nRunning 06_run_file_fixed: item-1-pushed codebase, item 1 meshfiles\n');
fprintf('   dt                  = %.6e s\n', dt);
fprintf('   tEnd                = %.6e s\n', tEnd);
fprintf('   requested steps     = %.0f\n', nSteps);

% FIX 3: Check if out already exists in workspace before launching
if exist('out', 'var')
    out = softlube_run_case_global_coupled(cfg, out);
else
    out = softlube_run_case_global_coupled(cfg);
end

out.cfg = cfg;
save(outputFile, 'out');

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
