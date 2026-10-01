function [cfg, gap1DWindow] = build_cfg_full2D_pressure2(varargin)
%BUILD_CFG_FULL2D_PRESSURE2  Same case as run_input_full2D_pressure2.m,
% refactored into a function so a grid-refinement study can call it
% repeatedly at different resolutions without duplicating/hand-editing
% the cfg block each time. All fields match run_input_full2D_pressure2.m
% exactly except for the resolution knobs, which are exposed as
% name-value arguments.
%
%   cfg = build_cfg_full2D_pressure2();
%   cfg = build_cfg_full2D_pressure2('NrExterior2D', 40, 'fineDz', 0.025e-6);
%
% Name-value overrides (defaults match run_input_full2D_pressure2.m):
%   NrExterior2D   (20)     radial cell count in each exterior 2D block
%   coarseDz       (1.0e-6) axial spacing away from the gap
%   fineDz         (0.05e-6) axial spacing inside the gap window
%   dt             (3e-4)   time step size [s]
%   nSteps         (1)      requested number of time steps

p = inputParser;
%added_temp
addParameter(p, 'useHybridGap1DExterior2DFluid', true);
%
addParameter(p, 'NrExterior2D', 20);
addParameter(p, 'coarseDz', 1.0e-6);
addParameter(p, 'fineDz', 0.05e-6);
addParameter(p, 'dt', 3e-4);
addParameter(p, 'nSteps', 1);
parse(p, varargin{:});
opt = p.Results;
%added_temp
useHybrid = opt.useHybridGap1DExterior2DFluid;
%
softlubeDir = '';

dt = opt.dt;
nSteps = opt.nSteps;
tEnd = nSteps * dt;
useAdaptiveTimeStep = true;
dtMin = dt / 1024;

zMin = -6e-6;
zMax = 10e-6;
rOuter = 4e-6;

pAtZMin = 0;
pAtZMax = 0;
axisPressureBC = 'dPdr=0';

gap1DWindow = [-0.2e-6, 4.2e-6];
finePressureWindow = gap1DWindow;
coarseDz = opt.coarseDz;
fineDz = opt.fineDz;

NrExterior2D = opt.NrExterior2D;
NrFluid2D = NrExterior2D;
radialFineWindow = [2e-6, 4e-6];
radialFineWeight = 4.0;

cfg = struct();
cfg.geometry.endotheliumPrestressFile = fullfile(softlubeDir, 'solid_endo_P300_wide.mat');
cfg.geometry.leukocytePrestressFile = fullfile(softlubeDir, 'solid_leu_P600.mat');
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
cfg.solid.leukocyte.deformable = true;
cfg.solid.leukocyte.EL = 200;
cfg.solid.leukocyte.nuL = 0.46;
cfg.solid.leukocyte.useViscoelasticLeukocyte = true;
cfg.solid.leukocyte.etaL = 1;
cfg.solid.leukocyte.supportL = 'axis';

cfg.solid.endothelium.Ee = 500;
cfg.solid.endothelium.nuE = 0.46;
cfg.solid.endothelium.useViscoelasticEndothelium = true;
cfg.solid.endothelium.etaE = 1;

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

% Plots/output are kept off here since this file is meant to be called
% many times in a loop by run_grid_refinement_study.m.
cfg.output.saveOutput = false;
cfg.output.outputFile = fullfile(softlubeDir, 'simulation_output_full2D_pressure.mat');
cfg.output.makePlots = false;
cfg.output.plotFinalGlobalPressureContour = false;
cfg.output.plotGlobal2DPressureTractionComparison = false;
cfg.output.printEvery = 1;
%added_temp
%cfg.output.storeFull2DFluidHist = false;
cfg.output.storeFull2DFluidHist = ~useHybrid;
%
cfg.output.store2DFluidArrays = true;
cfg.output.store2DPhysicalGridHist = true;
cfg.output.store2DSpeedHist = true;

cfg.runtime.useEnvironment = false;
cfg.ui.closeFigures = true;
cfg.parOverrides = struct();

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

cfg.fluid.useFull2DFluid = true;
%added_temp
%cfg.fluid.useHybridGap1DExterior2DFluid = true;
cfg.fluid.useHybridGap1DExterior2DFluid = useHybrid;
%
cfg.fluid.hybridGapZ = gap1DWindow;
cfg.fluid.hybridExteriorMinCells = 4;
cfg.fluid.hybridExteriorFailMode = 'warn';
cfg.fluid.hybridUseExterior2DPressureVector = true;
cfg.fluid.useBodyFittedMACFluid = true;
cfg.fluid.useBodyFittedMACTractionCorrection = true;
% Step 3 : correct until the worst interface mismatch
% is below a threshold, instead of a fixed pass count. See
% apply_bodyfitted_MAC_traction_correction_feedback.m.
cfg.fluid.useFeedbackTractionCorrection = true;
cfg.fluid.maxBodyFittedTractionCorrections = 1;
cfg.fluid.bodyFittedTractionCorrectionRelax = 0.05;
cfg.fluid.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.fluid.resolveFluidAfterTractionCorrection = true;
cfg.fluid.useBodyFittedMACTractionInSolid = true;

cfg.parOverrides.useFull2DFluid = true;
%added_temp
%cfg.parOverrides.useHybridGap1DExterior2DFluid = true;
cfg.parOverrides.useHybridGap1DExterior2DFluid = useHybrid;
%
cfg.parOverrides.hybridGapZ = gap1DWindow;
cfg.parOverrides.NrExterior2D = NrExterior2D;
cfg.parOverrides.bodyFittedRadialFineWindow = radialFineWindow;
cfg.parOverrides.bodyFittedRadialFineWeight = radialFineWeight;
cfg.parOverrides.hybridExteriorMinCells = 4;
cfg.parOverrides.hybridExteriorFailMode = 'warn';
cfg.parOverrides.hybridUseExterior2DPressureVector = true;
cfg.parOverrides.useBodyFittedMACFluid = true;
cfg.parOverrides.useBodyFittedMACTractionCorrection = true;
cfg.parOverrides.useFeedbackTractionCorrection = true;
cfg.parOverrides.maxBodyFittedTractionCorrections = 1;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = 0.05;
cfg.parOverrides.bodyFittedTractionCorrectionFailMode = 'warn';
cfg.parOverrides.resolveFluidAfterTractionCorrection = true;
cfg.parOverrides.useBodyFittedMACTractionInSolid = true;
cfg.parOverrides.NrFluid2D = NrFluid2D;
cfg.parOverrides.Nr = NrFluid2D;

cfg.parOverrides.useGlobal2DPressureTraction = false;

end
