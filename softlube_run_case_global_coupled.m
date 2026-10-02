
% =========================================================================
% SOFTLUBE_RUN_CASE_GLOBAL_COUPLED
% Coupled solver driver with global 2D body-fitted MAC fluid domain and
% axisymmetric finite-deformation solid mechanics.
% =========================================================================

function out = softlube_run_case_global_coupled(cfg,varargin)
%SOFTLUBE_RUN_CASE_GLOBAL_COUPLED Run coupled solver with a global pressure domain.
%   Includes numerical safeguards against element distortion and pressure runaway.

if nargin == 1
    if ~isfield(cfg, 'ui') || ~isfield(cfg.ui, 'closeFigures') || cfg.ui.closeFigures
        close all;
    end

    [par, S, uE_pre] = softlube_prepare_case(cfg);

    % STRICT GUARD: Full 2D Mode explicitly overrides hybrid/1D flags
    par.useFull2DFluid = true;
    par.useHybridGap1DExterior2DFluid = false;
    par.useGlobal2DPressureTraction = false;
    par.useSmoothHybridBlending = false;

    fprintf(['2D MAC/deformable-leukocyte mode: full2D=%s, ', ...
        'fixedCylinder=%s, prestressedLeukocyte=%s, exactInterface=%s, ', ...
        'global2DPressureTraction=%s, hybridGap1DExterior2D=%s.\n'], ...
        string_on_off(par.useFull2DFluid), ...
        string_on_off(par.useFixedCylindricalLeukocyte), ...
        string_on_off(par.usePrestressedLeukocyteIC), ...
        string_on_off(par.useExactDeformedInterface), ...
        string_on_off(use_global2d_pressure_traction(par)), ...
        string_on_off(use_hybrid_gap1d_exterior2d_fluid(par)));

    % Axial fluid grid (shared with interface interpolation locations)
    z = make_global_1d_z_grid(par);
    par.NzFluid = numel(z);
    par.zGrid = z;
    par.dz = mean(diff(z));

    % Building the finite-element meshes for the solid domains
    meshE = prepare_axisym_mesh_cache(S.meshE);
    interfaceE = S.interfaceE;
    baseE  = S.baseE;
    [deltaE_pre, ~,~] = exact_interface_radius_velocity( ...
        meshE, uE_pre, uE_pre, interfaceE, z, par.dt,par);

    % Leukocyte handling
    meshL = [];
    interfaceL = [];
    baseL = [];
    parL = [];
    uL_pre = [];
    deltaL_pre = par.RLout * ones(size(z));

    useFixedCylindricalLeukocyte = isfield(par,'useFixedCylindricalLeukocyte') && ...
        par.useFixedCylindricalLeukocyte;

    if ~(isfield(par, 'noLeukocyte') && par.noLeukocyte) && ~useFixedCylindricalLeukocyte
        SL = load(par.leukocytePrestressFile);
        preservedEL = par.EL;
        preservedNuL = par.nuL;
        par = apply_leukocyte_prestress_parameters(par, SL);
        par.EL = preservedEL;
        par.nuL = preservedNuL;
        par.GL = par.EL / (2 * (1 + par.nuL));
        par.KL = par.EL / (3 * (1 - 2 * par.nuL));
        meshL = SL.meshL;
        interfaceL = SL.interfaceL;
        baseL = SL.baseL;

        rigidLeukocyteAtInit = isfield(par, 'rigidLeukocyte') && par.rigidLeukocyte;
        useRLoutInnerAtInit = use_RLout_fluid_interface_for_solid_leukocyte(par);

        if ~par.usePrestressedLeukocyteIC
            if rigidLeukocyteAtInit && useRLoutInnerAtInit
                uL_pre = zeros(size(SL.uL_pre(:)));
                fprintf(['Rigid leukocyte: using fixed fluid inner boundary ', ...
                    'r = RLout, not prestressed leukocyte radius.\n']);
            else
                warning(['The loaded leukocyte prestress mesh is not valid with zero ', ...
                    'initial displacement for this coupled gap. Using uL_pre from %s.'], ...
                    par.leukocytePrestressFile);
                par.usePrestressedLeukocyteIC = true;
                uL_pre = SL.uL_pre(:);
                fprintf('Using prestressed leukocyte initial state from %s\n', ...
                    par.leukocytePrestressFile);
            end
        else
            uL_pre = SL.uL_pre(:);
            fprintf('Using prestressed leukocyte initial state from %s\n', ...
                par.leukocytePrestressFile);
        end

        if par.usePrestressedLeukocyteIC
            warn_leukocyte_prestress_load_mismatch(SL, par);
        end

        if useRLoutInnerAtInit
            deltaL_pre = par.RLout * ones(size(z));
        else
            [deltaL_pre, ~,~] = exact_interface_radius_velocity( ...
                meshL, uL_pre, uL_pre, interfaceL, z, par.dt,par);
            % Mask initial leukocyte boundary to prevent phantom overlap
            zL_nodes_init = meshL.nodes(:,2) + uL_pre(2:2:end);
            outOfLeukocyte_init = (z < min(zL_nodes_init)) | (z > max(zL_nodes_init));
            deltaL_pre(outOfLeukocyte_init) = par.RLout;
        end

        meshL = prepare_axisym_mesh_cache(meshL);
        parL = leukocyte_solid_parameters(par);
    else
        fprintf('Fixed cylindrical leukocyte: no leukocyte deformation, no leukocyte prestress; r_L = %.6e m.\n', par.RLout);
    end

    nStepsEnv = str2double(getenv('SOFTLUBE_NSTEPS'));
    if isfinite(nStepsEnv) && nStepsEnv > 0
        nSteps = max(1, round(nStepsEnv));
        par.tEnd = nSteps * par.dt;
    else
        nSteps = round(par.tEnd/par.dt);
    end

    historyCapacity = max(nSteps, 1);

    % Guarded History Vector Initialization
    tHist = nan(historyCapacity,1);
    dtHist = nan(historyCapacity,1);
    stepWallTimeHist = nan(historyCapacity,1);
    retryHist = zeros(historyCapacity,1);

    stateHist = cell(historyCapacity,1);
    fluidHist = cell(historyCapacity,1);
    deltaEHist = zeros(numel(z), historyCapacity);
    deltaLHist = zeros(numel(z), historyCapacity);
    pHist      = zeros(numel(z), historyCapacity);
    tauEHist   = zeros(numel(z), historyCapacity);
    tauLHist   = zeros(numel(z), historyCapacity);
    uzEHist    = zeros(numel(z), historyCapacity);
    uzLHist    = zeros(numel(z), historyCapacity);

    NrHist2D = par.NrFluid2D;
    if isfield(par,'Nr') && isfinite(par.Nr) && par.Nr > 0
        NrHist2D = par.Nr;
    end

    PHist      = nan(NrHist2D, numel(z), historyCapacity);
    urCHist    = nan(NrHist2D, numel(z), historyCapacity);
    uzCHist    = nan(NrHist2D, numel(z), historyCapacity);
    speedCHist = nan(NrHist2D, numel(z), historyCapacity);
    RPHist     = nan(NrHist2D, numel(z), historyCapacity);
    ZPHist     = nan(NrHist2D, numel(z), historyCapacity);
    p2DMaxHist = nan(historyCapacity,1);
    trEHist    = nan(2, numel(z), historyCapacity);
    trLHist    = nan(2, numel(z), historyCapacity);
    diagHist   = cell(historyCapacity,1);
    tractionCorrectionHistory = cell(historyCapacity,1);

    stoppedEarly = false;
    stopStep = 0;
    stopReason = '';

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        state = initial_state(z, par, meshE, uE_pre, deltaE_pre);
    else
        state = initial_state(z, par, meshE, uE_pre, deltaE_pre, ...
            meshL, interfaceL, uL_pre, deltaL_pre);
    end

    useExactDeformedInterface = isfield(par, 'useExactDeformedInterface') && ...
        par.useExactDeformedInterface;
    if useExactDeformedInterface
        if isfield(par, 'noLeukocyte') && par.noLeukocyte
            state = update_exact_fluid_interfaces( ...
                state, state, meshE, interfaceE, [], [], z, par);
        else
            state = update_exact_fluid_interfaces( ...
                state, state, meshE, interfaceE, meshL, interfaceL, z, par);
        end
    end
    if use_RLout_fluid_interface_for_solid_leukocyte(par)
        fprintf('Fluid inner boundary for leukocyte fixed at r = RLout = %.6e m.\n', par.RLout);
    end

    t0State = state;
    t0State.t = 0; % Explicit reference time anchor
    t0Fluid = [];
    try
        parT0 = par;
        if ~isfield(parT0, 'dt') || ~isfinite(parT0.dt) || parT0.dt <= 0
            parT0.dt = 1.0;
        end
        [fluid0, ok0, ~] = solve_selected_poststep_fluid(z, state, state, parT0);
        if ok0
            if isfield(fluid0, 'meshF') && ~isempty(fluid0.meshF) && ...
                    isfield(fluid0, 'ur2D') && ~isempty(fluid0.ur2D) && ...
                    isfield(fluid0, 'uz2D') && ~isempty(fluid0.uz2D)
                try
                    fluid0.meshF = add_fluid_nodes(fluid0.meshF);
                    [pCell0, sigmaCell0, center0] = recover_fluid_nodes_pressure_stress_Q4( ...
                        fluid0.meshF, fluid0.ur2D, fluid0.uz2D, par.mu, fluid0.pCell);
                    fluid0.pCellNode = pCell0;
                    fluid0.centerNode = center0;
                    fluid0.sigmaCellNode = sigmaCell0;
                catch
                    fluid0.pCellNode = [];
                    fluid0.centerNode = [];
                    fluid0.sigmaCellNode = [];
                end
            end
            t0Fluid = fluid0;
        else
            warning('t=0 reference fluid evaluation failed; out.t0Fluid will be empty.');
        end
    catch ME
        warning('t=0 reference fluid evaluation errored (%s); out.t0Fluid will be empty.', ME.message);
    end

    % Synchronize mass audit initial reference volume at t = 0
    z_vec0 = z(:);
    rE0 = t0State.deltaE(:);
    rL0 = t0State.deltaL(:);
    gap0 = rE0 - rL0;

    z_mid_target0 = 0.5 * (min(z_vec0) + max(z_vec0)); % Dynamic domain midplane
    [~, mid_idx0] = min(abs(z_vec0 - z_mid_target0));

    upper_mask0 = (1:numel(z_vec0))' > mid_idx0;
    upper_indices0 = find(upper_mask0);

    if ~isempty(upper_indices0)
        [~, min_upper_rel0] = min(gap0(upper_indices0));
        z_upper_min_idx0 = upper_indices0(min_upper_rel0);

        z_start_idx0 = mid_idx0;
        z_end_idx0 = z_upper_min_idx0;

        if z_end_idx0 > z_start_idx0
            dz_local0 = mean(diff(z_vec0));
            rE_cv0 = rE0(z_start_idx0:z_end_idx0);
            rL_cv0 = rL0(z_start_idx0:z_end_idx0);

            % Pre-calculate t = 0 reference upper-gap volume
            t0State.auditUpperVol = pi * sum((rE_cv0.^2 - rL_cv0.^2)) * dz_local0;
            state.auditUpperVol = t0State.auditUpperVol;
            state.t = 0;
        end
    end

    tn = 0;
    tNow = 0;
    dtNext = par.dt;
    timeTol = 100 * eps(max(par.tEnd, 1));

elseif nargin >= 2
    % =========================================================================
    % ARGUMENT 1: MAT FILE PATH OR LOADED STRUCT
    % ARGUMENT 2: TARGET TIMESTEP INDEX FOR ROLLBACK (e.g., 13)
    % =========================================================================
    fileOrData = cfg;
    targetStep = varargin{1};

    if nargin >= 3 && isstruct(varargin{2})
        runtimeCfg = varargin{2};
    elseif isstruct(cfg)
        runtimeCfg = cfg;
    else
        runtimeCfg = [];
    end

    % Load MAT file if a file path string was passed
    if ischar(fileOrData) || isstring(fileOrData)
        fprintf('\n[CHECKPOINT LOAD] Reading file: %s\n', fileOrData);
        S_check = load(fileOrData);
        if isfield(S_check, 'out')
            checkpointOut = S_check.out;
        else
            checkpointOut = S_check;
        end
    elseif isstruct(fileOrData)
        if isfield(fileOrData, 'out')
            checkpointOut = fileOrData.out;
        else
            checkpointOut = fileOrData;
        end
    else
        error('First input argument must be a MAT file path string or a loaded results struct.');
    end

    fprintf('\n=== [TRACE 1: MAT FILE CONTENTS] ===\n');
    if isfield(checkpointOut, 'par')
        fprintf('  - checkpointOut.par.useActiveTranslocation : %d\n', ...
            isfield(checkpointOut.par, 'useActiveTranslocation') && checkpointOut.par.useActiveTranslocation);
    end
    if isfield(checkpointOut, 'cfg') && isfield(checkpointOut.cfg, 'parOverrides')
        fprintf('  - checkpointOut.cfg.parOverrides.useActiveTranslocation : %d\n', ...
            isfield(checkpointOut.cfg.parOverrides, 'useActiveTranslocation') && checkpointOut.cfg.parOverrides.useActiveTranslocation);
    end

    % Cap target step to maximum available historical step in file
    maxAvailableStep = checkpointOut.stopStep;
    if isfield(checkpointOut, 'tHist') && ~isempty(checkpointOut.tHist)
        maxAvailableStep = nnz(isfinite(checkpointOut.tHist));
    end

    targetStep = min(round(targetStep), maxAvailableStep);

    fprintf('=== [CHECKPOINT RESUME] Extracting Step %d state (t = %.6e s) ===\n', ...
        targetStep, checkpointOut.t(targetStep));

    % Extract configuration and mesh structures from checkpoint
    par = checkpointOut.par;

    % =========================================================================
    % FORCE CHECKPOINT PARAMETERS ON RESTART (FIX 1)
    % =========================================================================
    if isstruct(cfg) && isfield(cfg, 'output') && isfield(cfg.output, 'checkpointFile') && ~isempty(cfg.output.checkpointFile)
        par.checkpointFile = cfg.output.checkpointFile;
        if isfield(cfg.output, 'checkpointEvery') && cfg.output.checkpointEvery > 0
            par.checkpointEvery = cfg.output.checkpointEvery;
        else
            par.checkpointEvery = 1;
        end
    elseif isfield(checkpointOut, 'cfg') && isfield(checkpointOut.cfg, 'output') && isfield(checkpointOut.cfg.output, 'checkpointFile')
        par.checkpointFile = checkpointOut.cfg.output.checkpointFile;
        par.checkpointEvery = 1;
    else
        % Default fallback if unassigned
        par.checkpointFile = fullfile(pwd, 'case_7_t.mat');
        par.checkpointEvery = 1;
    end

    % =========================================================================
    % [1002 FIX] RE-BASE STORED FILE REFERENCES ON THE CURRENT RUN
    % =========================================================================
    % par in the checkpoint carries absolute paths of the run that wrote it
    % (e.g. another user's or an older code folder). Write this run's output
    % next to its checkpoint in the current run folder, and read the
    % mesh/prestress files from the current code folder whenever the stored
    % path is not readable (otherwise baseL stays empty and the restart fails).
    % A checkpointFile/outputFile given in the runtime cfg still wins (below).
    [~, cpName1002, cpExt1002] = fileparts(par.checkpointFile);
    if isempty(cpName1002), cpName1002 = 'case_7_t'; cpExt1002 = '.mat'; end
    par.checkpointFile = fullfile(pwd, [cpName1002 cpExt1002]);
    par.outputFile = par.checkpointFile;
    codeDir1002 = fileparts(mfilename('fullpath'));
    pf1002 = {'leukocytePrestressFile', 'endotheliumPrestressFile'};
    for k1002 = 1:numel(pf1002)
        if isfield(par, pf1002{k1002}) && ~isempty(par.(pf1002{k1002})) && ...
                ~exist(par.(pf1002{k1002}), 'file')
            [~, nm1002, ex1002] = fileparts(par.(pf1002{k1002}));
            par.(pf1002{k1002}) = fullfile(codeDir1002, [nm1002 ex1002]);
        end
    end
    % [1002 FIX] Carry the t = 0 reference state through restarts (was dropped)
    if isfield(checkpointOut, 't0State'), t0State = checkpointOut.t0State; end
    if isfield(checkpointOut, 't0Fluid'), t0Fluid = checkpointOut.t0Fluid; end

    % =========================================================================
    % OVERRIDE LOADED PAR WITH RUNTIME CFG & PAROVERRIDES (RESTART FIX)
    % =========================================================================
    % 1. Apply parOverrides saved in checkpoint first
    % [1002 FIX] A finished run's .mat also contains out.cfg (added by the run
    % file). Its parOverrides are the run file's settings, not necessarily the
    % values the run used (out.par), and its file paths point to the old run
    % folder. Only fill in fields missing from par; never override the values
    % the run used and never take file locations from it. (The runtime cfg
    % passed to this call is still applied in full in step 2 below.)
    if isfield(checkpointOut, 'cfg') && isfield(checkpointOut.cfg, 'parOverrides')
        skip1002 = {'checkpointFile', 'outputFile', 'saveOutput', 'checkpointEvery'};
        fnames = fieldnames(checkpointOut.cfg.parOverrides);
        for k = 1:numel(fnames)
            if ~ismember(fnames{k}, skip1002) && ~isfield(par, fnames{k})
                par.(fnames{k}) = checkpointOut.cfg.parOverrides.(fnames{k});
            end
        end
    end

    % 2. Force fresh runtime cfg overrides passed in function call
    if isstruct(runtimeCfg)
        if isfield(runtimeCfg, 'parOverrides') && isstruct(runtimeCfg.parOverrides)
            fnames = fieldnames(runtimeCfg.parOverrides);
            for k = 1:numel(fnames)
                par.(fnames{k}) = runtimeCfg.parOverrides.(fnames{k});
            end
        end
        if isfield(runtimeCfg, 'numerics') && isstruct(runtimeCfg.numerics)
            fnames = fieldnames(runtimeCfg.numerics);
            for k = 1:numel(fnames)
                par.(fnames{k}) = runtimeCfg.numerics.(fnames{k});
            end
        end
        if isfield(runtimeCfg, 'solid') && isfield(runtimeCfg.solid, 'leukocyte') && isstruct(runtimeCfg.solid.leukocyte)
            fnames = fieldnames(runtimeCfg.solid.leukocyte);
            for k = 1:numel(fnames)
                par.(fnames{k}) = runtimeCfg.solid.leukocyte.(fnames{k});
            end
        end
    end
    % =========================================================================
    % =========================================================================


    if ~isfield(par, 'eta_solid') || isempty(par.eta_solid)
        par.eta_solid = 10.0; % [Pa*s] Transient solid viscosity default
    end

    % =========================================================================
    % STRATEGY A SAFEGUARD & STABILITY INITIALIZATION (GLOBAL DEFAULTS)
    % =========================================================================
    if ~isfield(par, 'maxFluidPressureCap') || isempty(par.maxFluidPressureCap)
        par.maxFluidPressureCap = 3000.0; % [Pa] Upper bound on fluid pressure
    end
    if ~isfield(par, 'maxFluidShearCap') || isempty(par.maxFluidShearCap)
        par.maxFluidShearCap = 500.0;     % [Pa] Upper bound on shear stress
    end
    if ~isfield(par, 'v_max_cap') || isempty(par.v_max_cap)
        par.v_max_cap = 0.020;            % [m/s] Kinematic velocity clamp (20 mm/s)
    end
    if ~isfield(par, 'solidAbsTol') || isempty(par.solidAbsTol)
        par.solidAbsTol = 1.0e-12;        % [N] Absolute force tolerance
    end
    if ~isfield(par, 'solidFallbackAbsTol') || isempty(par.solidFallbackAbsTol)
        par.solidFallbackAbsTol = 5.0e-9; % [N] Micro-element fallback floor (5 nN)
    end
    if ~isfield(par, 'newtonTolSolid') || isempty(par.newtonTolSolid)
        par.newtonTolSolid = 1.0e-3;      % Relative force tolerance
    end
    if ~isfield(par, 'eta_solid') || isempty(par.eta_solid)
        par.eta_solid = 10.0;             % [Pa*s] Transient solid viscosity
    end

    % Main Warm-Start Predictor Gain Factor
    alpha_pred = 0.30;
    % =========================================================================

    

    % =========================================================================
    % AUDIT PAR OVERRIDES TO IDENTIFY RESIDUAL SHIFT
    % =========================================================================
    fprintf('\n=== [PAR OVERRIDES AUDIT] ===\n');
    fprintf('  - par.dt          : %.6e s\n', par.dt);
    fprintf('  - par.uScaleMono  : %.6e\n', par.uScaleMono);
    fprintf('  - par.pScaleMono  : %.6e\n', par.pScaleMono);
    if isfield(par, 'eta_solid'), fprintf('  - par.eta_solid   : %.6e Pa*s\n', par.eta_solid); end
    if isfield(par, 'EL'),        fprintf('  - par.EL          : %.6e Pa\n', par.EL); end
    if isfield(par, 'nuL'),       fprintf('  - par.nuL         : %.6e\n', par.nuL); end
    fprintf('=============================\n\n');

    % =========================================================================
    % 2. EXTRACT STEP 13 STATE (t_n) & RESTORE CONVERGED FLUID FIELDS
    % =========================================================================
    if isfield(checkpointOut, 'stateHist') && numel(checkpointOut.stateHist) >= targetStep ...
            && ~isempty(checkpointOut.stateHist{targetStep})
        state = checkpointOut.stateHist{targetStep};
        stBaseTarget = checkpointOut.stateHist{targetStep};
    else
        state = checkpointOut.state;
        stBaseTarget = checkpointOut.state;
    end

    % DIRECT ASSIGNMENT: Force exact baseline uEPrev and uLPrev from checkpoint
    if isfield(stBaseTarget, 'uEPrev') && ~isempty(stBaseTarget.uEPrev)
        state.uEPrev = stBaseTarget.uEPrev;
    end
    if isfield(stBaseTarget, 'uLPrev') && ~isempty(stBaseTarget.uLPrev)
        state.uLPrev = stBaseTarget.uLPrev;
    end
    % =========================================================================
    % RESTORE 2D UNSTEADY STOKES FLUID FIELDS & TRACTIONS FROM FLUIDHIST
    % =========================================================================
    if isfield(checkpointOut, 'fluidHist') && numel(checkpointOut.fluidHist) >= targetStep ...
            && ~isempty(checkpointOut.fluidHist{targetStep})

        fl13 = checkpointOut.fluidHist{targetStep};

        % 1. 1D & 2D Tractions/Pressures (Directly synced to state)
        if isfield(fl13, 'p') && ~isempty(fl13.p)
            state.p = fl13.p(:);
            state.p2D = fl13.p(:);
        end
        if isfield(fl13, 'P') && ~isempty(fl13.P)
            state.P2DField = fl13.P;
        end
        if isfield(fl13, 'tauE') && ~isempty(fl13.tauE)
            state.tauE = fl13.tauE(:);
        end
        if isfield(fl13, 'tauL') && ~isempty(fl13.tauL)
            state.tauL = fl13.tauL(:);
        end

        % 2. 2D Face Velocities for Unsteady Momentum Term (du/dt)
        if isfield(fl13, 'ur') && ~isempty(fl13.ur)
            state.ur2DFaceField = fl13.ur;
        elseif isfield(fl13, 'ur2D') && ~isempty(fl13.ur2D)
            state.ur2DFaceField = fl13.ur2D;
        end

        if isfield(fl13, 'uz') && ~isempty(fl13.uz)
            state.uz2DFaceField = fl13.uz;
        elseif isfield(fl13, 'uz2D') && ~isempty(fl13.uz2D)
            state.uz2DFaceField = fl13.uz2D;
        end
    end

    % Extract or Reconstruct Step 13 (v^n) Velocities
    if isfield(state, 'v_wall_E') && ~isempty(state.v_wall_E)
        state.v_wall_E = state.v_wall_E;
    elseif isfield(checkpointOut, 'fluidHist') && numel(checkpointOut.fluidHist) >= targetStep ...
            && isfield(checkpointOut.fluidHist{targetStep}, 'v_wall_E')
        state.v_wall_E = checkpointOut.fluidHist{targetStep}.v_wall_E;
    elseif targetStep > 1 && isfield(checkpointOut, 'stateHist') && ~isempty(checkpointOut.stateHist{targetStep-1})
        [~, ~, state.v_wall_E] = exact_interface_radius_velocity( ...
            meshE, state.uE, checkpointOut.stateHist{targetStep-1}.uE, interfaceE, z, par.dt, par);
    end

    if isfield(state, 'v_wall_L') && ~isempty(state.v_wall_L)
        state.v_wall_L = state.v_wall_L;
    elseif isfield(checkpointOut, 'fluidHist') && numel(checkpointOut.fluidHist) >= targetStep ...
            && isfield(checkpointOut.fluidHist{targetStep}, 'v_wall_L')
        state.v_wall_L = checkpointOut.fluidHist{targetStep}.v_wall_L;
    end

    % Extract z-grid from checkpoint EARLY before initializing velocity buffers
    if isfield(checkpointOut, 'z') && ~isempty(checkpointOut.z)
        z = checkpointOut.z;
    else
        z = make_global_1d_z_grid(par);
    end
    % =========================================================================
    % 3. EXTRACT STEP 12 STATE (t_{n-1}) & RECONSTRUCT PREDICTOR BUFFERS (v^{n-1})
    % =========================================================================
    v_wall_E_nm1 = zeros(size(z));
    v_wall_L_nm1 = zeros(size(z));

    if targetStep > 1 && isfield(checkpointOut, 'stateHist') && ...
            numel(checkpointOut.stateHist) >= (targetStep - 1) && ...
            ~isempty(checkpointOut.stateHist{targetStep - 1})

        stateNm1 = checkpointOut.stateHist{targetStep - 1};

        % Viscoelastic displacement history for Step 14 (u^{13} - u^{12}) / dt
        % --- FIXED (EXACT BASELINE RECONSTRUCTION) ---
        % Reconstruct true uEPrev (Step 11 displacement) matching baseline stateHist{13}
        % --- GUARDED DISPLACEMENT HISTORY ---
        % Only fall back to stateNm2/stateNm1 if state.uEPrev or state.uLPrev was NOT
        % directly assigned from stBaseTarget (checkpointOut.stateHist{targetStep})
        if ~isfield(state, 'uEPrev') || isempty(state.uEPrev)
            if targetStep > 2 && isfield(checkpointOut, 'stateHist') && ...
                    numel(checkpointOut.stateHist) >= (targetStep - 2) && ...
                    ~isempty(checkpointOut.stateHist{targetStep - 2})
                stateNm2 = checkpointOut.stateHist{targetStep - 2};
                if isfield(stateNm2, 'uE') && ~isempty(stateNm2.uE), state.uEPrev = stateNm2.uE; end
            else
                if isfield(stateNm1, 'uE') && ~isempty(stateNm1.uE), state.uEPrev = stateNm1.uE; end
            end
        end

        if ~isfield(state, 'uLPrev') || isempty(state.uLPrev)
            if targetStep > 2 && isfield(checkpointOut, 'stateHist') && ...
                    numel(checkpointOut.stateHist) >= (targetStep - 2) && ...
                    ~isempty(checkpointOut.stateHist{targetStep - 2})
                stateNm2 = checkpointOut.stateHist{targetStep - 2};
                if isfield(stateNm2, 'uL') && ~isempty(stateNm2.uL), state.uLPrev = stateNm2.uL; end
            else
                if isfield(stateNm1, 'uL') && ~isempty(stateNm1.uL), state.uLPrev = stateNm1.uL; end
            end
        end
        if isfield(stateNm1, 'p')  && ~isempty(stateNm1.p),  state.pPrev  = stateNm1.p;  end

        % Step 12 (v^{n-1}) Endothelium Velocity Reconstruction
        if isfield(stateNm1, 'v_wall_E') && ~isempty(stateNm1.v_wall_E)
            v_wall_E_nm1 = stateNm1.v_wall_E;
        elseif isfield(checkpointOut, 'fluidHist') && numel(checkpointOut.fluidHist) >= (targetStep - 1) ...
                && isfield(checkpointOut.fluidHist{targetStep - 1}, 'v_wall_E')
            v_wall_E_nm1 = checkpointOut.fluidHist{targetStep - 1}.v_wall_E;
        elseif targetStep > 2 && isfield(checkpointOut.stateHist{targetStep-2}, 'uE')
            [~, ~, v_wall_E_nm1] = exact_interface_radius_velocity( ...
                meshE, stateNm1.uE, checkpointOut.stateHist{targetStep-2}.uE, interfaceE, z, par.dt, par);
        end

        % Step 12 (v^{n-1}) Leukocyte Velocity Reconstruction
        if isfield(stateNm1, 'v_wall_L') && ~isempty(stateNm1.v_wall_L)
            v_wall_L_nm1 = stateNm1.v_wall_L;
        elseif isfield(checkpointOut, 'fluidHist') && numel(checkpointOut.fluidHist) >= (targetStep - 1) ...
                && isfield(checkpointOut.fluidHist{targetStep - 1}, 'v_wall_L')
            v_wall_L_nm1 = checkpointOut.fluidHist{targetStep - 1}.v_wall_L;
        elseif targetStep > 2 && isfield(checkpointOut.stateHist{targetStep-2}, 'uL')
            if ~(isfield(par, 'noLeukocyte') && par.noLeukocyte) && ...
                    ~(isfield(par, 'useFixedCylindricalLeukocyte') && par.useFixedCylindricalLeukocyte) && ...
                    ~use_RLout_fluid_interface_for_solid_leukocyte(par)
                [~, ~, v_wall_L_nm1] = exact_interface_radius_velocity( ...
                    meshL, stateNm1.uL, checkpointOut.stateHist{targetStep-2}.uL, interfaceL, z, par.dt, par);
                zL_nodes = meshL.nodes(:,2) + stateNm1.uL(2:2:end);
                outOfLeukocyte = (z < min(zL_nodes)) | (z > max(zL_nodes));
                v_wall_L_nm1(outOfLeukocyte) = 0;
            end
        end
    end

    % Enable predictor extrapolation in solve_monolithic_two_solids_fsolve_timestep
    % [1002 FIX] Predictor state of the continuation: (1) dtPrev = dt actually used in the restart step
    % (the time loop stores stateNew.dtPrev = dt of the step; par.dt here made the fsolve warm start
    % y0 = u + (dt/dtPrev)(u - uPrev) different after a retried step); (2) a file written by the fixed
    % code (it contains baseE) keeps the predictor setting of the run that wrote it (a fresh start runs
    % without the predictor, so a restart from its file must too); older files (e.g. the step-13 input)
    % keep the original behaviour: predictor on.
    if isfield(checkpointOut, 'baseE') && isfield(checkpointOut, 'par') && isfield(checkpointOut.par, 'useMonoPredictor')
        par.useMonoPredictor = checkpointOut.par.useMonoPredictor;
    else
        par.useMonoPredictor = true;
    end
    state.dtPrev = par.dt;
    if isfield(checkpointOut, 'dtHist') && numel(checkpointOut.dtHist) >= targetStep && isfinite(checkpointOut.dtHist(targetStep))
        state.dtPrev = checkpointOut.dtHist(targetStep);
    end
    fprintf('[1002 FIX] Restart: useMonoPredictor = %d, dtPrev = %.4e s\n', par.useMonoPredictor, state.dtPrev);

    % =========================================================================
    % IN-LINE BASELINE VS RESTART PARITY CHECK
    % =========================================================================
    if targetStep <= numel(checkpointOut.stateHist) && ~isempty(checkpointOut.stateHist{targetStep})
        stBase = checkpointOut.stateHist{targetStep};
        fprintf('\n=== [BASELINE VS RESTART PARITY AUDIT] ===\n');
        if isfield(stBase, 'uE') && isfield(state, 'uE')
            fprintf('  1. ||uE_base - uE_restart|| : %.6e m\n', norm(stBase.uE - state.uE));
        end
        if isfield(stBase, 'uL') && isfield(state, 'uL')
            fprintf('  2. ||uL_base - uL_restart|| : %.6e m\n', norm(stBase.uL - state.uL));
        end
        if isfield(stBase, 'uEPrev') && isfield(state, 'uEPrev')
            fprintf('  3. ||uEPrev_base - uEPrev_restart|| : %.6e m\n', norm(stBase.uEPrev - state.uEPrev));
        end
        if isfield(stBase, 'uLPrev') && isfield(state, 'uLPrev')
            fprintf('  4. ||uLPrev_base - uLPrev_restart|| : %.6e m\n', norm(stBase.uLPrev - state.uLPrev));
        end
        if isfield(stBase, 'p') && isfield(state, 'p')
            fprintf('  5. ||p_base - p_restart||     : %.6e Pa\n', norm(stBase.p - state.p));
        end
        fprintf('===========================================\n\n');
    end

    nStepsEnv = str2double(getenv('SOFTLUBE_NSTEPS'));
    if isfinite(nStepsEnv) && nStepsEnv > 0
        nSteps = max(1, round(nStepsEnv));
        par.tEnd = nSteps * par.dt;
    else
        nSteps = round(par.tEnd/par.dt);
    end
    historyCapacity = max(nSteps, 1);

    % Reset time-stepping counters to roll back to Step 13 baseline
    tn = targetStep;
    tNow = state.t;
    dtNext = par.dt;
    % [1002 FIX] Continue with the time step the uninterrupted run would use next: the time loop sets
    % dtNext = min(par.dt, par.dtGrowFactor*dt) after a step that needed a retry or a reduced dt, but a
    % restart always started again from par.dt (e.g. step 15 at 2.5e-6 instead of 1.375e-6 after the
    % retried step 14), so the continuation differed from the uninterrupted run. Same rule, from the
    % saved history of the restart step.
    if isfield(par, 'enableAdaptiveTimeStep') && par.enableAdaptiveTimeStep && ...
            isfield(checkpointOut, 'dtHist') && numel(checkpointOut.dtHist) >= targetStep && ...
            isfield(checkpointOut, 'retryHist') && numel(checkpointOut.retryHist) >= targetStep
        dtPrev1002 = checkpointOut.dtHist(targetStep);
        if isfinite(dtPrev1002) && (checkpointOut.retryHist(targetStep) > 0 || dtPrev1002 < par.dt)
            dtNext = min(par.dt, par.dtGrowFactor * dtPrev1002);
            fprintf('[1002 FIX] Restart: next dt = %.4e s (step %d used dt = %.4e s with %d retries)\n', ...
                dtNext, targetStep, dtPrev1002, checkpointOut.retryHist(targetStep));
        end
    end
    timeTol = 100 * eps(max(par.tEnd, 1));
    stoppedEarly = false;
    stopStep = 0;
    stopReason = '';

    % Truncate history arrays cleanly up through targetStep
    tHist = checkpointOut.t(1:targetStep);
    tHist(targetStep+1:historyCapacity) = nan;

    dtHist = checkpointOut.dtHist(1:targetStep);
    dtHist(targetStep+1:historyCapacity) = nan;

    stepWallTimeHist = nan(historyCapacity,1);
    stepWallTimeHist(1:targetStep) = checkpointOut.stepWallTimeHist(1:targetStep);

    retryHist = zeros(historyCapacity,1);
    retryHist(1:targetStep) = checkpointOut.retryHist(1:targetStep);

    % =========================================================================
    % ROBUST BASE-E BOUNDARY NODE RESOLUTION
    % =========================================================================
    meshE = checkpointOut.meshE;
    interfaceE = checkpointOut.interfaceE;

    if isfield(checkpointOut, 'baseE') && ~isempty(checkpointOut.baseE)
        baseE = checkpointOut.baseE;
    elseif isfield(meshE, 'baseE') && ~isempty(meshE.baseE)
        baseE = meshE.baseE;
    else
        % Fallback: Load directly from endothelium prestress geometry file
        endoFile = '';
        if isfield(par, 'endotheliumPrestressFile') && exist(par.endotheliumPrestressFile, 'file')
            endoFile = par.endotheliumPrestressFile;
        elseif isstruct(cfg) && isfield(cfg, 'geometry') && isfield(cfg.geometry, 'endotheliumPrestressFile') && exist(cfg.geometry.endotheliumPrestressFile, 'file')
            endoFile = cfg.geometry.endotheliumPrestressFile;
        else
            endoFile = fullfile(fileparts(mfilename('fullpath')), 'solid_endothelium_P300.mat');
        end

        if exist(endoFile, 'file')
            fprintf('   [Resume Helper] Loading missing baseE boundary from %s\n', endoFile);
            SE = load(endoFile);
            baseE = SE.baseE;
        else
            error('softlube_run_case_global_coupled: baseE boundary nodes could not be resolved from checkpoint or %s', endoFile);
        end
    end

    % 0-based to 1-based Index Safeguard
    if any(baseE(:) == 0)
        baseE = baseE + 1;
        fprintf('   [Index Correction] Shifted baseE from 0-based to 1-based indexing.\n');
    end
    if exist('interfaceE', 'var') && ~isempty(interfaceE) && any(interfaceE(:) == 0)
        interfaceE = interfaceE + 1;
        fprintf('   [Index Correction] Shifted interfaceE from 0-based to 1-based indexing.\n');
    end
    meshE.baseE = baseE;

    hasLeukocyte = ~(isfield(par, 'noLeukocyte') && par.noLeukocyte) && ...
        ~(isfield(par, 'useFixedCylindricalLeukocyte') && par.useFixedCylindricalLeukocyte);

    meshL = []; interfaceL = []; baseL = []; parL = [];
    if hasLeukocyte
        if isfield(checkpointOut, 'meshL'), meshL = checkpointOut.meshL; end
        if isfield(checkpointOut, 'interfaceL'), interfaceL = checkpointOut.interfaceL; end
        if isfield(checkpointOut, 'baseL'), baseL = checkpointOut.baseL; end

        if isempty(baseL) && isfield(par, 'leukocytePrestressFile') && exist(par.leukocytePrestressFile, 'file')
            SL_load = load(par.leukocytePrestressFile);
            if isfield(SL_load, 'baseL'), baseL = SL_load.baseL; end
            if isfield(SL_load, 'interfaceL') && isempty(interfaceL), interfaceL = SL_load.interfaceL; end
            if isfield(SL_load, 'meshL') && isempty(meshL), meshL = SL_load.meshL; end
        end

        if ~isempty(baseL) && any(baseL(:) == 0)
            baseL = baseL + 1;
            fprintf('   [Index Correction] Shifted baseL from 0-based to 1-based indexing.\n');
        end
        if ~isempty(interfaceL) && any(interfaceL(:) == 0)
            interfaceL = interfaceL + 1;
            fprintf('   [Index Correction] Shifted interfaceL from 0-based to 1-based indexing.\n');
        end
        if ~isempty(baseL) && isstruct(meshL)
            meshL.baseL = baseL;
        end

        parL = leukocyte_solid_parameters(par);
    end

    z = checkpointOut.z;

    % Truncate cell history arrays up to targetStep
    stateHist = cell(historyCapacity,1);
    stateHist(1:targetStep) = checkpointOut.stateHist(1:targetStep);

    fluidHist = cell(historyCapacity,1);
    fluidHist(1:targetStep) = checkpointOut.fluidHist(1:targetStep);

    deltaEHist = zeros(numel(z), historyCapacity);
    deltaEHist(:, 1:targetStep) = checkpointOut.deltaEHist(:, 1:targetStep);

    deltaLHist = zeros(numel(z), historyCapacity);
    deltaLHist(:, 1:targetStep) = checkpointOut.deltaLHist(:, 1:targetStep);

    pHist = zeros(numel(z), historyCapacity);
    pHist(:, 1:targetStep) = checkpointOut.pHist(:, 1:targetStep);

    tauEHist = zeros(numel(z), historyCapacity);
    tauEHist(:, 1:targetStep) = checkpointOut.tauEHist(:, 1:targetStep);

    tauLHist = zeros(numel(z), historyCapacity);
    tauLHist(:, 1:targetStep) = checkpointOut.tauLHist(:, 1:targetStep);

    uzEHist = zeros(numel(z), historyCapacity);
    uzEHist(:, 1:targetStep) = checkpointOut.uzEHist(:, 1:targetStep);

    uzLHist = zeros(numel(z), historyCapacity);
    uzLHist(:, 1:targetStep) = checkpointOut.uzLHist(:, 1:targetStep);

    NrHist2D = par.NrFluid2D;
    if isfield(par, 'Nr') && isfinite(par.Nr) && par.Nr > 0, NrHist2D = par.Nr; end

    PHist      = nan(NrHist2D, numel(z), historyCapacity);
    urCHist    = nan(NrHist2D, numel(z), historyCapacity);
    uzCHist    = nan(NrHist2D, numel(z), historyCapacity);
    speedCHist = nan(NrHist2D, numel(z), historyCapacity);
    RPHist     = nan(NrHist2D, numel(z), historyCapacity);
    ZPHist     = nan(NrHist2D, numel(z), historyCapacity);

    if isfield(checkpointOut, 'PHist') && size(checkpointOut.PHist, 3) >= targetStep
        PHist(:,:,1:targetStep) = checkpointOut.PHist(:,:,1:targetStep);
    end
    if isfield(checkpointOut, 'urCHist') && size(checkpointOut.urCHist, 3) >= targetStep
        urCHist(:,:,1:targetStep) = checkpointOut.urCHist(:,:,1:targetStep);
    end
    if isfield(checkpointOut, 'uzCHist') && size(checkpointOut.uzCHist, 3) >= targetStep
        uzCHist(:,:,1:targetStep) = checkpointOut.uzCHist(:,:,1:targetStep);
    end
    if isfield(checkpointOut, 'speedCHist') && size(checkpointOut.speedCHist, 3) >= targetStep
        speedCHist(:,:,1:targetStep) = checkpointOut.speedCHist(:,:,1:targetStep);
    end
    if isfield(checkpointOut, 'RPHist') && size(checkpointOut.RPHist, 3) >= targetStep
        RPHist(:,:,1:targetStep) = checkpointOut.RPHist(:,:,1:targetStep);
    end
    if isfield(checkpointOut, 'ZPHist') && size(checkpointOut.ZPHist, 3) >= targetStep
        ZPHist(:,:,1:targetStep) = checkpointOut.ZPHist(:,:,1:targetStep);
    end

    p2DMaxHist = nan(historyCapacity,1);
    p2DMaxHist(1:targetStep) = checkpointOut.p2DMaxHist(1:targetStep);

    trEHist = nan(2, numel(z), historyCapacity);
    trEHist(:,:,1:targetStep) = checkpointOut.trEHist(:,:,1:targetStep);

    trLHist = nan(2, numel(z), historyCapacity);
    trLHist(:,:,1:targetStep) = checkpointOut.trLHist(:,:,1:targetStep);

    diagHist = cell(historyCapacity,1);
    diagHist(1:targetStep) = checkpointOut.diagHist(1:targetStep);

    tractionCorrectionHistory = cell(historyCapacity,1);
    if isfield(checkpointOut, 'tractionCorrectionHistory') && numel(checkpointOut.tractionCorrectionHistory) >= targetStep
        tractionCorrectionHistory(1:targetStep) = checkpointOut.tractionCorrectionHistory(1:targetStep);
    end
end
% =========================================================================
% GUARANTEE TIME-STEPPING SCOPE VARIABLES FOR BOTH FRESH & RESTART MODES
% =========================================================================

if ~exist('par', 'var') || isempty(par)
    if exist('checkpointOut', 'var') && isfield(checkpointOut, 'par')
        par = checkpointOut.par;
    else
        error('softlube_run_case_global_coupled: Structure "par" could not be resolved.');
    end
end

if ~exist('tNow', 'var') || isempty(tNow)
    if exist('state', 'var') && isfield(state, 't') && ~isempty(state.t)
        tNow = state.t;
    else
        tNow = 0.0;
    end
end

if ~exist('tn', 'var') || isempty(tn)
    if exist('targetStep', 'var') && ~isempty(targetStep)
        tn = targetStep;
    else
        tn = 0;
    end
end

if ~exist('timeTol', 'var') || isempty(timeTol)
    timeTol = 100 * eps(max(par.tEnd, 1));
end

if ~exist('dtNext', 'var') || isempty(dtNext)
    dtNext = par.dt;
end
% =========================================================================


% Main Time-Stepping Loop
while tNow < par.tEnd - timeTol
    old = state;
    stepTicId = tic;

    rigidLeukocyte = isfield(par, 'rigidLeukocyte') && par.rigidLeukocyte;
    hasLeukocyte = ~(isfield(par, 'noLeukocyte') && par.noLeukocyte);

    dtAttempt = min(dtNext, par.tEnd - tNow);
    retryCount = 0;
    acceptedStep = false;
    stepReason = '';

    % =========================================================================
    % [CHECKPOINT 1] STRATEGY A: DAMPED KINEMATIC WARM-START PREDICTOR
    % =========================================================================
    if tn == 0 || ~isfield(old, 'v_wall_E') || isempty(old.v_wall_E)
        v_wall_E_step_init = zeros(size(z));
        v_wall_L_step_init = zeros(size(z));
    elseif tn == 1
        v_wall_E_step_init = old.v_wall_E + alpha_pred * (old.v_wall_E - zeros(size(z)));
        v_wall_L_step_init = old.v_wall_L + alpha_pred * (old.v_wall_L - zeros(size(z)));
    else
        v_wall_E_step_init = old.v_wall_E + alpha_pred * (old.v_wall_E - v_wall_E_nm1);
        v_wall_L_step_init = old.v_wall_L + alpha_pred * (old.v_wall_L - v_wall_L_nm1);
    end

    if tn > 0
        norm_v_prev = norm(old.v_wall_E);
        norm_v_pred = norm(v_wall_E_step_init);
        fprintf('  [Strategy A Predictor] Step %d Warm-Start: ||v_n|| = %.3e | ||v_pred|| = %.3e (alpha = %.2f)\n', ...
            tn + 1, norm_v_prev, norm_v_pred, alpha_pred);
    end

    fprintf('\n[Step Init Audit] Step %d Baseline Anchored: max|v_E_init| = %.3e m/s, max|v_L_init| = %.3e m/s\n', ...
        tn + 1, max(abs(v_wall_E_step_init)), max(abs(v_wall_L_step_init)));

    while ~acceptedStep

        parStep = par;
        parStep.dt = dtAttempt;

        v_wall_E_last = v_wall_E_step_init;
        v_wall_L_last = v_wall_L_step_init;

        parStep.noLeukocyte                  = false;
        parStep.useFixedCylindricalLeukocyte = false;
        parStep.rigidLeukocyte               = false;
        if isstruct(parL)
            parL.dt = dtAttempt;
        end

        parStep.useFull2DFluid = true;
        parStep.useHybridGap1DExterior2DFluid = false;
        parStep.useGlobal2DPressureTraction = false;
        parStep.useSmoothHybridBlending = false;

        parStep.useGapRepulsion = false;
        parStep.gapFloor = 1.0e-8;

        if isfield(parStep, 'fluid2DPenaltyFactor')
            parStep.penaltyLambda = parStep.fluid2DPenaltyFactor * ...
                parStep.mu / max(parStep.dt, realmin);
        end

        tNew = tNow + dtAttempt;

        try
            if isfield(old, 'geometryE') && isfield(old.geometryE, 'JEmin')
                if old.geometryE.JEmin < 0.25
                    fprintf('  [Mesh Safeguard] Endothelium minJ = %.4e < 0.25; applying nodal relaxation.\n', old.geometryE.JEmin);
                    meshE = relax_surface_mesh_nodes(meshE, 0.10);
                end
            end

            maxCouplingIters = 35;
            couplingTol      = 1e-2;
            omega_base       = 0.25;

            omega_min  = 0.08;
            omega_max  = 0.50;
            omega_init = 0.08;
            alpha_pred = 0.30;

            if ~exist('v_wall_E_nm1', 'var') || isempty(v_wall_E_nm1)
                v_wall_E_nm1 = zeros(size(z));
            end
            if ~exist('v_wall_L_nm1', 'var') || isempty(v_wall_L_nm1)
                v_wall_L_nm1 = zeros(size(z));
            end

            dV_max = 1.0e-3;
            omega_effective = 0.05;
            omega_old       = omega_effective;
            R_prev          = [];

            stateIterInput = old;

            if isfield(stateIterInput, 'p2D') && ~isempty(stateIterInput.p2D) && isvector(stateIterInput.p2D)
                T_norm_applied = stateIterInput.p2D(:);
            else
                T_norm_applied = stateIterInput.p(:);
            end

            if isfield(stateIterInput, 'tauE') && ~isempty(stateIterInput.tauE)
                T_tang_applied = stateIterInput.tauE(:);
            else
                T_tang_applied = zeros(size(T_norm_applied));
            end

            pRelaxedFinal = old.p;

            for couplingIter = 1:maxCouplingIters

                if couplingIter == 1
                    stateIterInput.uE = old.uE;
                    stateIterInput.p  = old.p(:);
                    if isfield(old, 'uEPrev') && ~isempty(old.uEPrev)
                        stateIterInput.uEPrev = old.uEPrev;
                    else
                        stateIterInput.uEPrev = old.uE;
                    end

                    if isfield(old, 'uL')
                        stateIterInput.uL = old.uL;
                        if isfield(old, 'uLPrev') && ~isempty(old.uLPrev)
                            stateIterInput.uLPrev = old.uLPrev;
                        else
                            stateIterInput.uLPrev = old.uL;
                        end
                    end
                end

                if hasLeukocyte && ~rigidLeukocyte
                    stateTrial = solve_monolithic_two_solids_fsolve_timestep( ...
                        stateIterInput, meshE, interfaceE, baseE, ...
                        meshL, interfaceL, baseL, parL, ...
                        z, parStep);
                else
                    stateTrial = solve_monolithic_analytical_timestep( ...
                        stateIterInput, meshE, interfaceE, baseE, z, parStep);
                end

                useExactDeformedInterface = isfield(parStep, 'useExactDeformedInterface') && ...
                    parStep.useExactDeformedInterface;

                if useExactDeformedInterface
                    if isfield(parStep, 'noLeukocyte') && parStep.noLeukocyte
                        stateTrial = update_exact_fluid_interfaces( ...
                            stateTrial, old, meshE, interfaceE, [], [], z, parStep);
                    else
                        stateTrial = update_exact_fluid_interfaces( ...
                            stateTrial, old, meshE, interfaceE, meshL, interfaceL, z, parStep);
                    end
                end

                [stateTrial.deltaE, stateTrial.UwE, v_E_raw] = exact_interface_radius_velocity( ...
                    meshE, stateTrial.uE, old.uE, interfaceE, z, parStep.dt, parStep);

                v_L_raw = zeros(size(z));
                if hasLeukocyte && ~isempty(meshL) && ~(isfield(parStep, 'noLeukocyte') && parStep.noLeukocyte)
                    if use_RLout_fluid_interface_for_solid_leukocyte(parStep)
                        stateTrial.deltaL = parStep.RLout * ones(size(z));
                        stateTrial.UwL = zeros(size(z));
                        v_L_raw = zeros(size(z));
                    else
                        [stateTrial.deltaL, stateTrial.UwL, v_L_raw] = exact_interface_radius_velocity( ...
                            meshL, stateTrial.uL, old.uL, interfaceL, z, parStep.dt, parStep);

                        zL_nodes_trial = meshL.nodes(:,2) + stateTrial.uL(2:2:end);
                        outOfLeukocyte_trial = (z < min(zL_nodes_trial)) | (z > max(zL_nodes_trial));
                        stateTrial.deltaL(outOfLeukocyte_trial) = parStep.RLout;
                        stateTrial.UwL(outOfLeukocyte_trial) = 0;
                        v_L_raw(outOfLeukocyte_trial) = 0;
                    end
                end

                if hasLeukocyte && ~isempty(meshL)
                    stateTrial = attach_state_geometry_checks(stateTrial, meshE, meshL, parStep);
                else
                    stateTrial = attach_state_geometry_checks(stateTrial, meshE, [], parStep);
                end

                stateTrial.dt = parStep.dt;

                diff_v_E = v_E_raw - v_wall_E_last;
                delta_v_E_capped = sign(diff_v_E) .* min(abs(diff_v_E), dV_max);

                stateTrial.v_wall_E = v_wall_E_last + delta_v_E_capped;
                parStep.v_wall_E    = stateTrial.v_wall_E;

                if hasLeukocyte && ~isempty(meshL) && ~(isfield(parStep, 'noLeukocyte') && parStep.noLeukocyte) ...
                        && ~use_RLout_fluid_interface_for_solid_leukocyte(parStep)
                    diff_v_L = v_L_raw - v_wall_L_last;
                    delta_v_L_capped = sign(diff_v_L) .* min(abs(diff_v_L), dV_max);

                    stateTrial.v_wall_L = v_wall_L_last + delta_v_L_capped;
                    parStep.v_wall_L    = stateTrial.v_wall_L;
                else
                    stateTrial.v_wall_L = zeros(size(z));
                    parStep.v_wall_L    = zeros(size(z));
                end

                fprintf('    -> [Iter %02d Kinematic Audit] Endothelium: max|v_raw| = %.3e m/s | max|v_clamped| = %.3e m/s | max|dv_k| = %.3e m/s\n', ...
                    couplingIter, max(abs(v_E_raw)), max(abs(stateTrial.v_wall_E)), max(abs(diff_v_E)));

                v_wall_E_last = stateTrial.v_wall_E;
                v_wall_L_last = stateTrial.v_wall_L;

                if couplingIter == 1
                    fprintf('\n=== [STEP 14 ITER 1 INPUT AUDIT] ===\n');
                    fprintf('  - min(deltaE)  : %.8e m\n', min(stateTrial.deltaE));
                    fprintf('  - min(deltaL)  : %.8e m\n', min(stateTrial.deltaL));
                    fprintf('  - min(gap)     : %.8e m\n', min(stateTrial.deltaE - stateTrial.deltaL));
                    fprintf('  - max|v_wall_E|: %.8e m/s\n', max(abs(parStep.v_wall_E)));
                    fprintf('  - max|v_wall_L|: %.8e m/s\n', max(abs(parStep.v_wall_L)));
                    fprintf('  - z-grid points: %d (dz = %.3e m)\n', numel(z), parStep.dz);
                    fprintf('=====================================\n\n');
                end

                % =========================================================================
                % GAPFLOOR GUARDED FALLBACK CHECK
                % Safeguard against gap collapse (h < gapFloor) before calling 2D fluid solver
                % =========================================================================
                if isfield(parStep, 'gapFloor') && ~isempty(parStep.gapFloor)
                    gap_floor = parStep.gapFloor;
                else
                    gap_floor = 1.0e-8; % Default 10 nm floor fallback
                end

                gapTrial = stateTrial.deltaE(:) - stateTrial.deltaL(:);
                minGapFound = min(gapTrial);

                if minGapFound < gap_floor
                    fprintf('    [Gap Safeguard] Minimum gap (%.4e m) below floor (%.4e m). Enforcing floor.\n', ...
                        minGapFound, gap_floor);

                    % Clamp endothelium interface to preserve physical gap floor
                    stateTrial.deltaE = max(stateTrial.deltaE(:), stateTrial.deltaL(:) + gap_floor);
                end

                [fluidTrial, okFluid, fluidReason] = ...
                    solve_selected_poststep_fluid(z, old, stateTrial, parStep);
                if ~okFluid
                    error('Fluid solve failed at iteration %d: %s', couplingIter, fluidReason);
                end

                max_normal_traction_cap = parStep.maxFluidPressureCap;

                % =========================================================================
                % C1 CONTINUOUS (TANH) SMOOTH FLUID TRACTION SATURATION
                % =========================================================================
                pCap = parStep.maxFluidPressureCap; % 3000.0 Pa default
                if isfield(fluidTrial, 'p') && ~isempty(fluidTrial.p)
                    max_p_raw = max(abs(fluidTrial.p(:)));
                    if max_p_raw > pCap
                        fprintf('    [Traction Guard - C1 Smooth] Peak pressure (%.2f Pa > %.2f Pa). Applying smooth tanh saturation.\n', ...
                            max_p_raw, pCap);

                        % Smooth saturation function: p_sat = pCap * tanh(p / pCap)
                        fluidTrial.p = pCap * tanh(fluidTrial.p / pCap);
                        if isfield(fluidTrial, 'P') && ~isempty(fluidTrial.P)
                            fluidTrial.P = pCap * tanh(fluidTrial.P / pCap);
                        end
                    end
                end

                tauCap = parStep.maxFluidShearCap; % 500.0 Pa default
                if isfield(fluidTrial, 'tauE') && ~isempty(fluidTrial.tauE)
                    if max(abs(fluidTrial.tauE(:))) > tauCap
                        fluidTrial.tauE = tauCap * tanh(fluidTrial.tauE / tauCap);
                    end
                end
                if isfield(fluidTrial, 'tauL') && ~isempty(fluidTrial.tauL)
                    if max(abs(fluidTrial.tauL(:))) > tauCap
                        fluidTrial.tauL = tauCap * tanh(fluidTrial.tauL / tauCap);
                    end
                end
                % =========================================================================

                if isfield(parStep,'useFull2DFluid') && parStep.useFull2DFluid && ...
                        isfield(parStep,'useBodyFittedMACTractionCorrection') && ...
                        parStep.useBodyFittedMACTractionCorrection && ...
                        fluid_supports_partitioned_traction_correction(fluidTrial)
                    meshLcorr = []; interfaceLcorr = []; baseLcorr = []; parLcorr = [];
                    if hasLeukocyte && exist('meshL','var') && exist('interfaceL','var') && ...
                            exist('baseL','var') && exist('parL','var')
                        meshLcorr = meshL; interfaceLcorr = interfaceL; baseLcorr = baseL; parLcorr = parL;
                    end

                    [stateTrial, fluidTrial, okFluid, fluidReason] = ...
                        apply_bodyfitted_MAC_traction_correction( ...
                        z, old, stateTrial, fluidTrial, ...
                        meshE, interfaceE, baseE, ...
                        meshLcorr, interfaceLcorr, baseLcorr, parLcorr, ...
                        parStep);

                    if ~okFluid
                        error('Body-fitted MAC traction correction failed at iteration %d: %s', couplingIter, fluidReason);
                    end
                end

                if isfield(fluidTrial, 'p') && ~isempty(fluidTrial.p)
                    T_norm_fluid = fluidTrial.p(:);
                else
                    T_norm_fluid = zeros(size(z(:)));
                end

                if isfield(fluidTrial, 'tauE') && ~isempty(fluidTrial.tauE)
                    T_tang_fluid = fluidTrial.tauE(:);
                else
                    T_tang_fluid = zeros(size(T_norm_fluid));
                end

                nNodesCheck = min([numel(T_norm_fluid), numel(T_norm_applied), numel(T_tang_fluid), numel(T_tang_applied)]);

                if nNodesCheck > 0
                    mismatch_norm = T_norm_fluid(1:nNodesCheck) - T_norm_applied(1:nNodesCheck);
                    mismatch_tang = T_tang_fluid(1:nNodesCheck) - T_tang_applied(1:nNodesCheck);

                    L2_stress_mismatch  = sqrt(mean(mismatch_norm.^2 + mismatch_tang.^2));
                    Max_stress_mismatch = max(sqrt(mismatch_norm.^2 + mismatch_tang.^2));
                else
                    L2_stress_mismatch  = NaN;
                    Max_stress_mismatch = NaN;
                end

                relDiff = max(abs(T_norm_fluid - T_norm_applied)) / max(max(abs(T_norm_fluid)), 1.0);

                R_curr = T_norm_fluid - T_norm_applied;

                if couplingIter == 1
                    omega_raw       = omega_init;
                    omega_effective = omega_init;
                elseif couplingIter == 2 || isempty(R_prev)
                    omega_raw       = omega_init;
                    omega_effective = omega_init;
                else
                    delta_R = R_curr - R_prev;
                    denom = sum(delta_R.^2);

                    if denom > 1e-20
                        mu_k = - omega_old * (sum(R_prev .* delta_R) / denom);
                        omega_raw = omega_old + mu_k;

                        omega_filtered = 0.30 * omega_raw + 0.70 * omega_old;
                        omega_effective = max(omega_min, min(omega_max, omega_filtered));
                    else
                        omega_raw       = omega_min;
                        omega_effective = omega_min;
                    end

                    if omega_raw < parStep.omega_min
                        fprintf('    [Aitken Controller] Iter %02d: Raw omega (%.4f) hit FLOOR clamp -> Enforced omega = %.4f\n', ...
                            couplingIter, omega_raw, omega_effective);
                    elseif omega_raw > parStep.omega_max
                        fprintf('    [Aitken Controller] Iter %02d: Raw omega (%.4f) hit CEILING clamp -> Enforced omega = %.4f\n', ...
                            couplingIter, omega_raw, omega_effective);
                    else
                        fprintf('    [Aitken Controller] Iter %02d: Dynamic update active -> omega = %.4f (relDiff = %.3e)\n', ...
                            couplingIter, omega_effective, relDiff);
                    end
                end

                R_prev    = R_curr;
                omega_old = omega_effective;

                if omega_raw < omega_min
                    fprintf('    [Strategy A Aitken] Iter %02d: Raw omega (%.4f) hit FLOOR clamp -> Enforced omega = %.4f\n', ...
                        couplingIter, omega_raw, omega_effective);
                elseif omega_raw > omega_max
                    fprintf('    [Strategy A Aitken] Iter %02d: Raw omega (%.4f) hit CEILING clamp -> Enforced omega = %.4f\n', ...
                        couplingIter, omega_raw, omega_effective);
                else
                    fprintf('    [Strategy A Aitken] Iter %02d: Dynamic update active -> omega = %.4f\n', ...
                        couplingIter, omega_effective);
                end

                fprintf('    -> Coupling Iter %d/%d: relDiff = %.4e (tol = %.1e, omega_Aitken = %.4f)\n', ...
                    couplingIter, maxCouplingIters, relDiff, couplingTol, omega_effective);
                fprintf('       [FSI Verification] L2 Mismatch = %.4e Pa | Max Mismatch = %.4e Pa\n', ...
                    L2_stress_mismatch, Max_stress_mismatch);
                fprintf('       [FSI Verification] Max Fluid Norm = %.4e Pa | Max Applied Norm = %.4e Pa\n', ...
                    max(abs(T_norm_fluid)), max(abs(T_norm_applied)));

                if couplingIter == 1 || mod(couplingIter, 5) == 0 || relDiff < couplingTol
                    maxP    = max(abs(T_norm_fluid));
                    maxTauE = max(abs(T_tang_fluid));
                    fprintf('    [Intermediate Summary Iter %02d] Max(P) = %.2f Pa | Max(tauE) = %.2f Pa | relDiff = %.4e\n', ...
                        couplingIter, maxP, maxTauE, relDiff);
                end

                if relDiff < couplingTol
                    fprintf('  [Traction Coupling] Converged at iteration %d (relDiff = %.3e)\n', couplingIter, relDiff);
                    pRelaxedFinal = T_norm_fluid;
                    break;
                end

                % Automatic step-halving fallback if Aitken coupling diverges severely
                if couplingIter > 10 && relDiff > 1.5
                    error('Aitken:Divergence', ...
                        'Coupling divergence detected (relDiff = %.3e > 1.5). Forcing time-step reduction.', relDiff);
                end

                T_norm_applied = omega_effective * T_norm_fluid + (1.0 - omega_effective) * T_norm_applied;
                T_tang_applied = omega_effective * T_tang_fluid + (1.0 - omega_effective) * T_tang_applied;
                pRelaxedFinal  = T_norm_applied;

                stateIterInput = old;
                stateIterInput.p        = T_norm_applied;
                stateIterInput.p2D      = T_norm_applied;
                stateIterInput.pReduced = T_norm_applied;

                stateIterInput.tauE     = T_tang_applied;
                if isfield(fluidTrial, 'tauL') && ~isempty(fluidTrial.tauL)
                    stateIterInput.tauL = fluidTrial.tauL(:);
                else
                    stateIterInput.tauL = T_tang_applied;
                end

                stateIterInput.uE = stateTrial.uE;
                stateIterInput.uEPrev = old.uE;
                if isfield(stateTrial, 'uL')
                    stateIterInput.uL = stateTrial.uL;
                    stateIterInput.uLPrev = old.uL;
                end

                if isfield(fluidTrial, 'P') && ~isempty(fluidTrial.P)
                    stateIterInput.P2DField = fluidTrial.P;
                end
            end

            stateTrial.v_wall_E = parStep.v_wall_E;
            stateTrial.v_wall_L = parStep.v_wall_L;

            fprintf('\n  [Post-Coupling Handoff Audit] Step %d (t=%.4e s):\n', tn + 1, tNow + dtAttempt);
            fprintf('    -> Endothelium (v_wall_E) : max|v| = %.3e m/s, mean = %+.3e m/s\n', ...
                max(abs(stateTrial.v_wall_E)), mean(stateTrial.v_wall_E));
            if hasLeukocyte && ~(isfield(parStep, 'noLeukocyte') && parStep.noLeukocyte) && ~use_RLout_fluid_interface_for_solid_leukocyte(parStep)
                fprintf('    -> Leukocyte   (v_wall_L) : max|v| = %.3e m/s, mean = %+.3e m/s\n\n', ...
                    max(abs(stateTrial.v_wall_L)), mean(stateTrial.v_wall_L));
            else
                fprintf('    -> Leukocyte   (v_wall_L) : FIXED / BYPASSED (v_wall_L = 0)\n\n');
            end

            if isfield(parStep,'useFull2DFluid') && parStep.useFull2DFluid
                stateTrial.pReduced = stateTrial.p;
                stateTrial.p = pRelaxedFinal;
                stateTrial.p2D = pRelaxedFinal;
                stateTrial.pL2D = fluidTrial.pL;
                stateTrial.pE2D = fluidTrial.pE;
                if isfield(fluidTrial,'P') && ~isempty(fluidTrial.P)
                    stateTrial.P2DField = fluidTrial.P;
                end
                if isfield(fluidTrial,'urC') && ~isempty(fluidTrial.urC)
                    stateTrial.urC2DField = fluidTrial.urC;
                end
                if isfield(fluidTrial,'uzC') && ~isempty(fluidTrial.uzC)
                    stateTrial.uzC2DField = fluidTrial.uzC;
                end
                if isfield(fluidTrial,'ur') && ~isempty(fluidTrial.ur)
                    stateTrial.ur2DFaceField = fluidTrial.ur;
                end
                if isfield(fluidTrial,'uz') && ~isempty(fluidTrial.uz)
                    stateTrial.uz2DFaceField = fluidTrial.uz;
                end
                if isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes
                    if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Rur')
                        stateTrial.Rur2DField = fluidTrial.meshF.Rur;
                    end
                    if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Ruz')
                        stateTrial.Ruz2DField = fluidTrial.meshF.Ruz;
                    end
                end
                if isfield(fluidTrial,'meshF') && isfield(fluidTrial.meshF,'Rp')
                    stateTrial.Rp2DField = fluidTrial.meshF.Rp;
                    stateTrial.Zp2DField = fluidTrial.meshF.Zp;
                end
                if isfield(fluidTrial,'tractionE'), stateTrial.tractionE2D = fluidTrial.tractionE; end
                if isfield(fluidTrial,'tractionL'), stateTrial.tractionL2D = fluidTrial.tractionL; end
            else
                stateTrial.p = pRelaxedFinal;
                stateTrial.pReduced = pRelaxedFinal;
            end
            stateTrial.tauE = fluidTrial.tauE;
            stateTrial.tauL = fluidTrial.tauL;

            [stateTrial, fluidTrial, pressureLimited, pressureLimitReason] = ...
                apply_pressure_temporal_limiter(stateTrial, fluidTrial, old, z, parStep);

            if exist('fluidTrial', 'var') && isfield(fluidTrial, 'p')
                p_field = fluidTrial.p;
                if isfield(fluidTrial, 'P') && ~isempty(fluidTrial.P)
                    p_field = fluidTrial.P;
                end

                negative_mask = p_field < 0;
                if any(negative_mask(:))
                    neg_count = nnz(negative_mask);
                    total_count = numel(p_field);
                    neg_fraction = (neg_count / total_count) * 100;
                    min_neg_val = min(p_field(negative_mask));
                    max_pos_val = max(p_field(~negative_mask));

                    fprintf(['[Negative Pressure Audit] Step %d (t=%.4e s): ' ...
                        'Found negative pressure in %.2f%% of fluid domain ' ...
                        '(Min P = %.3e Pa, Max P = %.3e Pa).\n'], ...
                        tn + 1, tNow + dtAttempt, neg_fraction, min_neg_val, max_pos_val);

                    if neg_fraction > 15.0
                        fprintf('[Cavitation Safeguard] Truncating unphysical negative pressures to zero absolute pressure.\n');
                        p_field(negative_mask) = 0;
                        if isfield(fluidTrial, 'P') && ~isempty(fluidTrial.P)
                            fluidTrial.P = p_field;
                        else
                            fluidTrial.p = p_field;
                        end
                    end
                end
            end

            if exist('fluidTrial', 'var') && isfield(stateTrial, 'deltaE') && isfield(stateTrial, 'deltaL')
                z_vec = z(:);
                rE = stateTrial.deltaE(:);
                rL = stateTrial.deltaL(:);
                gap = rE - rL;

                z_mid_target = 0.5 * (min(z_vec) + max(z_vec));
                [~, mid_idx] = min(abs(z_vec - z_mid_target));

                upper_mask = (1:numel(z_vec))' > mid_idx;
                upper_indices = find(upper_mask);

                if ~isempty(upper_indices)
                    [~, min_upper_rel] = min(gap(upper_indices));
                    z_start_idx = mid_idx;
                    z_end_idx = upper_indices(min_upper_rel);

                    if z_end_idx > z_start_idx
                        dz_local = mean(diff(z_vec));
                        dz_cv = dz_local;
                        rE_cv = rE(z_start_idx:z_end_idx);
                        rL_cv = rL(z_start_idx:z_end_idx);

                        current_upper_volume = pi * sum((rE_cv.^2 - rL_cv.^2)) * dz_local;

                        if isfield(old, 'auditUpperVol') && isfield(old, 't')
                            prev_upper_vol_val = old.auditUpperVol;
                            prev_audit_t_val = old.t;
                        else
                            prev_upper_vol_val = current_upper_volume;
                            prev_audit_t_val = tNow;
                        end

                        stateTrial.auditUpperVol = current_upper_volume;

                        dt_step = (tNow + dtAttempt) - prev_audit_t_val;

                        if dt_step > 0
                            d_vol_upper_dt = (current_upper_volume - prev_upper_vol_val) / dt_step;

                            if isfield(fluidTrial, 'Q') && ~isempty(fluidTrial.Q)
                                Q_faces = fluidTrial.Q(:);

                                j_mid   = z_start_idx;
                                j_upper = min(z_end_idx + 1, numel(Q_faces));

                                Q_mid   = Q_faces(j_mid);
                                Q_upper = Q_faces(j_upper);

                                net_open_flux = Q_upper-Q_mid;

                                rE_sub = rE(z_start_idx:z_end_idx);
                                wrE_sub = parStep.v_wall_E(z_start_idx:z_end_idx);
                                Q_kin_E = sum(2 * pi * rE_sub .* wrE_sub) * dz_cv;

                                Q_kin_L = 0;
                                if hasLeukocyte && ~(isfield(parStep, 'noLeukocyte') && parStep.noLeukocyte)&& ~use_RLout_fluid_interface_for_solid_leukocyte(parStep)
                                    rL_sub = rL(z_start_idx:z_end_idx);
                                    wrL_sub = parStep.v_wall_L(z_start_idx:z_end_idx);
                                    Q_kin_L = -sum(2 * pi * rL_sub .* wrL_sub) * dz_cv;
                                end

                                Q_kinematic = Q_kin_E + Q_kin_L;

                                upper_mass_residual = abs(d_vol_upper_dt - (Q_kinematic - net_open_flux));
                                scale_Q = max([abs(d_vol_upper_dt), abs(Q_kinematic), abs(net_open_flux), realmin]);
                                rel_residual = upper_mass_residual / scale_Q;

                                if (upper_mass_residual > 1e-14 || rel_residual > 0.01) && dt_step > 1e-7
                                    fprintf(['[Upper Mass Audit] Step %d (t=%.4e s): ' ...
                                        'dV_upper/dt = %.3e m^3/s, Q_kin = %.3e m^3/s, Q_open = %.3e m^3/s | ' ...
                                        'Res = %.3e m^3/s, RelRes = %.2f%%\n'], ...
                                        tn + 1, tNow + dtAttempt, d_vol_upper_dt, Q_kinematic, net_open_flux, ...
                                        upper_mass_residual, rel_residual * 100);
                                end
                            end
                        end
                    end
                end
            end

            check_pressure_jump_retry(old, fluidTrial, z, parStep, tn);

            acceptedStep = true;
        catch ME
            stepReason = ['Monolithic solve failed: ', ME.message];

            if should_retry_time_step(stepReason, dtAttempt, retryCount, par)
                dtNew = max(par.dtMin, par.dtRetryFactor * dtAttempt);
                fprintf(['   retrying time step from t=%.6e s: ', ...
                    'dt %.3e -> %.3e after %s\n'], ...
                    tNow, dtAttempt, dtNew, compact_failure_reason(stepReason));
                dtAttempt = dtNew;
                retryCount = retryCount + 1;
                continue;
            end

            stoppedEarly = true;
            stopStep = tn;
            stopReason = stepReason;
            warning('Stopped early after step %d, attempted t = %.6e s. %s', ...
                tn, tNow + dtAttempt, stopReason);
            break;
        end
    end

    if stoppedEarly
        break;
    end

    tn = tn + 1;
    if tn > historyCapacity
        growBy = max(historyCapacity, max(nSteps, 1));
        newCapacity = historyCapacity + growBy;
        tHist(newCapacity,1) = nan;
        dtHist(newCapacity,1) = nan;
        stepWallTimeHist(newCapacity,1) = nan;
        retryHist(newCapacity,1) = 0;
        stateHist{newCapacity,1} = [];
        fluidHist{newCapacity,1} = [];
        deltaEHist(:,newCapacity) = 0;
        deltaLHist(:,newCapacity) = 0;
        pHist(:,newCapacity) = 0;
        tauEHist(:,newCapacity) = 0;
        tauLHist(:,newCapacity) = 0;
        uzEHist(:,newCapacity) = 0;
        uzLHist(:,newCapacity) = 0;

        PHist = cat(3, PHist, nan(NrHist2D, numel(z), growBy));
        urCHist = cat(3, urCHist, nan(NrHist2D, numel(z), growBy));
        uzCHist = cat(3, uzCHist, nan(NrHist2D, numel(z), growBy));
        speedCHist = cat(3, speedCHist, nan(NrHist2D, numel(z), growBy));
        RPHist = cat(3, RPHist, nan(NrHist2D, numel(z), growBy));
        ZPHist = cat(3, ZPHist, nan(NrHist2D, numel(z), growBy));

        p2DMaxHist(newCapacity,1) = nan;
        trEHist = cat(3, trEHist, nan(2, numel(z), growBy));
        trLHist = cat(3, trLHist, nan(2, numel(z), growBy));
        diagHist{newCapacity,1} = [];
        tractionCorrectionHistory{newCapacity,1} = [];
        historyCapacity = newCapacity;
    end

    if tn > 1
        v_wall_E_nm1 = old.v_wall_E;
        v_wall_L_nm1 = old.v_wall_L;
        fprintf('  [Strategy A Buffer] Accepted Step %d: Advanced v_wall_E_nm1 buffer (norm = %.3e)\n', ...
            tn + 1, norm(v_wall_E_nm1));
    end

    state = stateTrial;
    fluid = fluidTrial;
    tNow = tNew;
    state.t = tNow;

    maxSolidChange = max(abs(state.uE - old.uE));
    if isfield(state, 'uL') && ~isempty(state.uL)
        maxSolidChange = max(maxSolidChange, max(abs(state.uL - old.uL)));
    end
    oldFluidPressure = old.p;
    if isfield(old, 'p2D') && isvector(old.p2D)
        oldFluidPressure = old.p2D;
    end
    maxPressureChange = max(abs(fluid.p(:) - oldFluidPressure(:)));
    printEvery = 1;
    if isfield(par, 'printEvery') && isfinite(par.printEvery) && par.printEvery > 0
        printEvery = par.printEvery;
    end
    doPrintStep = tn == 1 || tNow >= par.tEnd - timeTol || ...
        mod(tn, printEvery) == 0 || retryCount > 0;
    if doPrintStep

        fprintf(['Time step %d (t=%.4e, dt=%.3e): min(deltaE)=%.6e, ', ...
            'max|p|=%.6e, max|du|=%.3e, max|dp|=%.3e, iters=%d, retries=%d\n'], ...
            tn, tNow, dtAttempt, min(state.deltaE), max(abs(fluid.p)), ...
            maxSolidChange, maxPressureChange, couplingIter, retryCount);

    end

    fluidStore = fluid;
    if isfield(par,'storeFull2DFluidHist') && ~par.storeFull2DFluidHist && ...
            isfield(fluidStore,'meshType') && strcmpi(fluidStore.meshType,'bodyfitted_MAC')
        fluidStore.ur = [];
        fluidStore.uz = [];
        if isfield(fluidStore,'meshF') && isfield(fluidStore.meshF,'Rp')
            meshLight = struct();
            meshLight.Rp = fluidStore.meshF.Rp;
            meshLight.Zp = fluidStore.meshF.Zp;
            meshLight.zc = fluidStore.meshF.zc;
            meshLight.zF = fluidStore.meshF.zF;
            if isfield(fluidStore.meshF,'Nr') && isfield(fluidStore.meshF,'Nz')
                meshLight.Nr = fluidStore.meshF.Nr;
                meshLight.Nz = fluidStore.meshF.Nz;
            else
                [meshLight.Nr, meshLight.Nz] = size(fluidStore.meshF.Rp);
            end
            if isfield(fluidStore.meshF,'deltaL_c')
                meshLight.deltaL_c = fluidStore.meshF.deltaL_c;
            end
            if isfield(fluidStore.meshF,'deltaE_c')
                meshLight.deltaE_c = fluidStore.meshF.deltaE_c;
            end
            fluidStore.meshF = meshLight;
        end
    end

    if isfield(fluidStore, 'meshF') && ~isempty(fluidStore.meshF) && ...
            isfield(fluidStore, 'ur2D') && ~isempty(fluidStore.ur2D) && ...
            isfield(fluidStore, 'uz2D') && ~isempty(fluidStore.uz2D)
        try
            fluidStore.meshF = add_fluid_nodes(fluidStore.meshF);
            [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4( ...
                fluidStore.meshF, fluidStore.ur2D, fluidStore.uz2D, par.mu, fluidStore.pCell);

            fluidStore.pCellNode = pCell;
            fluidStore.centerNode = center;
            fluidStore.sigmaCellNode = sigmaCell;
        catch
            fluidStore.pCellNode = [];
            fluidStore.centerNode = [];
            fluidStore.sigmaCellNode = [];
        end
    end

    stateHist{tn} = state;
    fluidHist{tn} = fluidStore;
    tHist(tn) = tNow;
    dtHist(tn) = dtAttempt;
    stepWallTimeHist(tn) = toc(stepTicId);
    retryHist(tn) = retryCount;
    deltaEHist(:,tn) = state.deltaE;
    deltaLHist(:,tn) = state.deltaL;
    pHist(:,tn)      = fluid.p;
    tauEHist(:,tn)   = fluid.tauE;
    tauLHist(:,tn)   = fluid.tauL;
    uzEHist(:,tn)    = fluid.uzE;
    uzLHist(:,tn)    = fluid.uzL;

    if isfield(fluid,'P') && ~isempty(fluid.P)
        [nrP,nzP] = size(fluid.P);
        nrS = min(size(PHist,1), nrP);
        nzS = min(size(PHist,2), nzP);
        PHist(1:nrS,1:nzS,tn) = fluid.P(1:nrS,1:nzS);
    end
    if isfield(fluid,'urC') && ~isempty(fluid.urC)
        [nrU,nzU] = size(fluid.urC);
        nrS = min(size(urCHist,1), nrU);
        nzS = min(size(urCHist,2), nzU);
        urCHist(1:nrS,1:nzS,tn) = fluid.urC(1:nrS,1:nzS);
    end
    if isfield(fluid,'uzC') && ~isempty(fluid.uzC)
        [nrU,nzU] = size(fluid.uzC);
        nrS = min(size(uzCHist,1), nrU);
        nzS = min(size(uzCHist,2), nzU);
        uzCHist(1:nrS,1:nzS,tn) = fluid.uzC(1:nrS,1:nzS);
    end
    if isfield(fluid,'urC') && isfield(fluid,'uzC') && ...
            ~isempty(fluid.urC) && ~isempty(fluid.uzC)
        [nrU,nzU] = size(fluid.urC);
        nrS = min(size(speedCHist,1), nrU);
        nzS = min(size(speedCHist,2), nzU);
        speedCHist(1:nrS,1:nzS,tn) = sqrt(fluid.urC(1:nrS,1:nzS).^2 + fluid.uzC(1:nrS,1:nzS).^2);
    end
    if isfield(fluid,'P') && ~isempty(fluid.P)
        Pvals = abs(fluid.P(:));
        Pvals = Pvals(isfinite(Pvals));
        if ~isempty(Pvals), p2DMaxHist(tn) = max(Pvals); end
    end
    if isfield(fluid,'meshF') && isfield(fluid.meshF,'Rp') && ~isempty(fluid.meshF.Rp)
        [nrR,nzR] = size(fluid.meshF.Rp);
        nrS = min(size(RPHist,1), nrR);
        nzS = min(size(RPHist,2), nzR);
        RPHist(1:nrS,1:nzS,tn) = fluid.meshF.Rp(1:nrR,1:nzR);
    end
    if isfield(fluid,'meshF') && isfield(fluid.meshF,'Zp') && ~isempty(fluid.meshF.Zp)
        [nrZ,nzZ] = size(fluid.meshF.Zp);
        nrS = min(size(ZPHist,1), nrZ);
        nzS = min(size(ZPHist,2), nzZ);
        ZPHist(1:nrZ,1:nzZ,tn) = fluid.meshF.Zp(1:nrZ,1:nzZ);
    end
    if isfield(fluid,'tractionE')
        trEHist(1,:,tn) = safe_interp1_same_or_resample(fluid.tractionE.z(:), fluid.tractionE.normal(:), z(:), 'fluid.tractionE.normal');
        trEHist(2,:,tn) = safe_interp1_same_or_resample(fluid.tractionE.z(:), fluid.tractionE.tangent(:), z(:), 'fluid.tractionE.tangent');
    end
    if isfield(fluid,'tractionL')
        trLHist(1,:,tn) = safe_interp1_same_or_resample(fluid.tractionL.z(:), fluid.tractionL.normal(:), z(:), 'fluid.tractionL.normal');
        trLHist(2,:,tn) = safe_interp1_same_or_resample(fluid.tractionL.z(:), fluid.tractionL.tangent(:), z(:), 'fluid.tractionL.tangent');
    end

    diag = compute_step_diagnostics(z, old, state, fluid, parStep, ...
        pressureLimited, pressureLimitReason);
    diag.dt = dtAttempt;
    diag.retries = retryCount;
    diagHist{tn} = diag;

    if isfield(par, 'checkpointFile') && ~isempty(par.checkpointFile) && ...
            isfield(par, 'checkpointEvery') && par.checkpointEvery > 0 && ...
            mod(tn, par.checkpointEvery) == 0 && tn > 0
        try
            out = struct();
            out.z = z;
            out.t = tHist(isfinite(tHist));
            out.dtHist = dtHist(isfinite(dtHist));
            out.stepWallTimeHist = stepWallTimeHist(1:tn);
            out.retryHist = retryHist(1:tn);
            if exist('t0State', 'var')
                out.t0State = t0State;
                out.t0Fluid = t0Fluid;
            end
            out.state = state;
            out.stateHist = stateHist(1:tn);
            out.fluidHist = fluidHist(1:tn);
            out.deltaEHist = deltaEHist(:,1:tn);
            out.deltaLHist = deltaLHist(:,1:tn);
            out.pHist = pHist(:,1:tn);
            out.tauEHist = tauEHist(:,1:tn);
            out.tauLHist = tauLHist(:,1:tn);
            out.uzEHist = uzEHist(:,1:tn);
            out.uzLHist = uzLHist(:,1:tn);
            out.PHist = PHist(:,:,1:tn);
            out.urCHist = urCHist(:,:,1:tn);
            out.uzCHist = uzCHist(:,:,1:tn);
            out.speedCHist = speedCHist(:,:,1:tn);
            out.RPHist = RPHist(:,:,1:tn);
            out.ZPHist = ZPHist(:,:,1:tn);
            out.p2DMaxHist = p2DMaxHist(1:tn);
            out.trEHist = trEHist(:,:,1:tn);
            out.trLHist = trLHist(:,:,1:tn);
            out.diagHist = diagHist(1:tn);
            out.tractionCorrectionHistory = tractionCorrectionHistory(1:tn);
            out.par = par;
            out.baseE = baseE;                                   % [1002 FIX]
            if exist('baseL','var') && ~isempty(baseL), out.baseL = baseL; end  % [1002 FIX]
            out.meshE = meshE;
            out.interfaceE = interfaceE;
            if exist('meshL','var') && ~isempty(meshL)
                out.meshL = meshL;
            end
            if exist('interfaceL','var') && ~isempty(interfaceL)
                out.interfaceL = interfaceL;
            end
            out.stoppedEarly = false;
            out.stopStep = tn;
            out.stopReason = '';
            out.isCheckpoint = true;
            cpDir = fileparts(par.checkpointFile);
            if ~isempty(cpDir) && ~exist(cpDir, 'dir')
                mkdir(cpDir);
            end
            save(par.checkpointFile, 'out', '-v7.3');
            clear out;
        catch MEcp
            warning('Checkpoint save failed at step %d (continuing run): %s', tn, MEcp.message);
        end
    end

    gapMinNow = diag.gapMin;
    if doPrintStep
        fprintf('   min gap after accepted step = %.6e m\n', gapMinNow);
        if isfield(par, 'diagnosticsEnabled') && par.diagnosticsEnabled
            print_step_diagnostics(diag);
        end
    end

    warn_step_diagnostics(diag, par);

    if isfield(par, 'stopAtMinGap') && par.stopAtMinGap && ...
            gapMinNow <= par.gapStopFactor * parStep.minGap
        stoppedEarly = true;
        stopStep = tn;
        stopReason = sprintf('Reached minimum gap: hmin = %.6e m', gapMinNow);
        fprintf('Reached minimum gap at step %d, t = %.6e s. hmin = %.6e m\n', ...
            tn, tNew, gapMinNow);
        break;
    end

    if isfield(par, 'enableAdaptiveTimeStep') && par.enableAdaptiveTimeStep && ...
            (retryCount > 0 || dtAttempt < par.dt)
        dtNext = min(par.dt, par.dtGrowFactor * dtAttempt);
    else
        dtNext = par.dt;
    end

    % [1002] Optional short runs: stop after a given step WITHOUT changing
    % par.tEnd (SOFTLUBE_NSTEPS shortens tEnd, which makes the last step use
    % dt = tEnd - t and stores the short tEnd in the .mat file).
    % Only active when SOFTLUBE_STOP_AFTER_STEP is set.
    stopAfter1002 = str2double(getenv('SOFTLUBE_STOP_AFTER_STEP'));
    if isfinite(stopAfter1002) && tn >= stopAfter1002
        fprintf('[STOP AFTER STEP] Stopping after step %d (SOFTLUBE_STOP_AFTER_STEP).\n', tn);
        break;
    end
end

if ~stoppedEarly
    stopStep = tn;
end

if stopStep < 1
    safeStep = 1;
else
    safeStep = stopStep;
end

stateHist = stateHist(1:stopStep);
fluidHist = fluidHist(1:stopStep);
deltaEHist = deltaEHist(:,1:stopStep);
deltaLHist = deltaLHist(:,1:stopStep);
pHist      = pHist(:,1:stopStep);
tauEHist   = tauEHist(:,1:stopStep);
tauLHist   = tauLHist(:,1:stopStep);
uzEHist    = uzEHist(:,1:stopStep);
uzLHist    = uzLHist(:,1:stopStep);
PHist      = PHist(:,:,1:stopStep);
urCHist    = urCHist(:,:,1:stopStep);
uzCHist    = uzCHist(:,:,1:stopStep);
speedCHist = speedCHist(:,:,1:stopStep);
RPHist     = RPHist(:,:,1:stopStep);
ZPHist     = ZPHist(:,:,1:stopStep);
p2DMaxHist = p2DMaxHist(1:stopStep);
trEHist    = trEHist(:,:,1:stopStep);
trLHist    = trLHist(:,:,1:stopStep);
diagHist   = diagHist(1:stopStep);
tractionCorrectionHistory = tractionCorrectionHistory(1:stopStep);

validSteps = find(isfinite(tHist));
if ~isempty(validSteps)
    lastValid = validSteps(end);
else
    lastValid = safeStep;
end

tHist            = tHist(1:lastValid);
dtHist           = dtHist(1:lastValid);
stepWallTimeHist = stepWallTimeHist(1:stopStep);
retryHist        = retryHist(1:stopStep);

adaptiveSummary = summarize_time_step_adaptation( ...
    dtHist, retryHist, par, stoppedEarly, stopReason);
if adaptiveSummary.totalRetries > 0 || stoppedEarly
    print_adaptive_summary(adaptiveSummary);
end

out = struct();
out.z = z;
out.t = tHist;
out.dtHist = dtHist;
out.stepWallTimeHist = stepWallTimeHist;
out.retryHist = retryHist;
if exist('t0State', 'var')
    out.t0State = t0State;
    out.t0Fluid = t0Fluid;
end
out.adaptiveSummary = adaptiveSummary;
out.state = state;
out.stateHist = stateHist;
out.fluidHist = fluidHist;
out.deltaEHist = deltaEHist;
out.deltaLHist = deltaLHist;
out.pHist = pHist;
out.tauEHist = tauEHist;
out.tauLHist = tauLHist;
out.uzEHist = uzEHist;
out.uzLHist = uzLHist;
out.PHist = PHist;
out.urCHist = urCHist;
out.uzCHist = uzCHist;
out.speedCHist = speedCHist;
out.RPHist = RPHist;
out.ZPHist = ZPHist;
out.p2DMaxHist = p2DMaxHist;
out.trEHist = trEHist;
out.trLHist = trLHist;
out.diagHist = diagHist;
out.tractionCorrectionHistory = tractionCorrectionHistory;
out.par = par;
out.baseE = baseE;                                   % [1002 FIX]
if exist('baseL','var') && ~isempty(baseL), out.baseL = baseL; end  % [1002 FIX]
out.meshE = meshE;
out.interfaceE = interfaceE;
if exist('meshL','var') && ~isempty(meshL)
    out.meshL = meshL;
end
if exist('interfaceL','var') && ~isempty(interfaceL)
    out.interfaceL = interfaceL;
end

out.stoppedEarly = stoppedEarly;
out.stopStep = stopStep;
out.stopReason = stopReason;

if isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure
    out.global1D = build_global_1d_pressure_view(out.z, out.pHist, par);
    if out.stopStep >= 1
        out.global1D = add_global_1d_blank_solid_pressure_view(out, par);
    end
end

if out.stopStep >= 1 && use_global2d_pressure_traction(par)
    out.global2DPressureTractionComparison = ...
        build_global2d_pressure_traction_comparison(out);
end

if out.stopStep < 1
    warning('Simulation failed before completing the first time step. No plots generated.');
    return;
end

nPlot = out.stopStep;
statePlot = out.stateHist{nPlot};
fluidPlot = out.fluidHist{nPlot};

if par.saveOutput
    outputFile = 'simulation_output_two_solid_full_analytical_nopre.mat';
    if isfield(par, 'outputFile') && ~isempty(par.outputFile)
        outputFile = par.outputFile;
    end
    save(outputFile,'out','-v7.3');
end

plotFinalGlobalPressureContour = isfield(par, 'useGlobal1DPressure') && ...
    par.useGlobal1DPressure && ...
    (par.makePlots || (isfield(par, 'plotFinalGlobalPressureContour') && ...
    par.plotFinalGlobalPressureContour));
if plotFinalGlobalPressureContour
    plot_global_1d_pressure_contour_with_blank_solids(out, par);
end

plotGlobal2DPressureTractionComparison = use_global2d_pressure_traction(par) && ...
    isfield(out, 'global2DPressureTractionComparison') && ...
    out.global2DPressureTractionComparison.available && ...
    (par.makePlots || (isfield(par, 'plotGlobal2DPressureTractionComparison') && ...
    par.plotGlobal2DPressureTractionComparison));
if plotGlobal2DPressureTractionComparison
    plot_global2d_pressure_traction_comparison(out, par);
end

if ~par.makePlots
    return;
end

[R, Z, Uz] = velocity_field_for_plot(z, fluidPlot, statePlot, par);

figure;
set(gca, 'FontSize', 24);
plot(z*1e6, fluidPlot.p, 'LineWidth', 1.8);
grid off;
xlabel('z [\mum]');
ylabel('Pressure [Pa]');
title(sprintf('Pressure field at t = %.4f s', out.t(nPlot)));

figure;
zMidIdx = round(numel(out.z) / 2);
plot(out.t, out.pHist(zMidIdx, :), 'LineWidth', 2);
set(gca, 'FontSize', 24);
xlabel('t [s]');
ylabel('Pressure at mid-point [Pa]');
title(sprintf('Pressure evolution at z = %.3f \\mum', out.z(zMidIdx)*1e6));
grid off;

if isfield(out,'PHist') && ~isempty(out.PHist)
    Pplot = out.PHist(:,:,nPlot);
    if any(isfinite(Pplot(:)))
        figure;
        set(gca, 'FontSize', 24);
        contourf(out.RPHist(:,:,nPlot)*1e6, out.ZPHist(:,:,nPlot)*1e6, ...
            Pplot, 40, 'LineColor', 'none');
        colorbar;
        hold on;
        plot(statePlot.deltaL*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        plot(statePlot.deltaE*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
        xlabel('r [\mum]');
        ylabel('z [\mum]');
        title(sprintf('Native 2D pressure P(r,z) at t = %.4f s', out.t(nPlot)));
    end
end

end

% ========================================================================
% HELPER FUNCTIONS
% ========================================================================
function tf = fluid_supports_partitioned_traction_correction(fluid)
tf = isstruct(fluid) && ...
    isfield(fluid, 'meshType') && ...
    strcmpi(fluid.meshType, 'bodyfitted_MAC') && ...
    isfield(fluid, 'p') && ~isempty(fluid.p);
end

function g1D = build_global_1d_pressure_view(z, pHist, par)
g1D = struct();
g1D.z = z(:);
g1D.pHist = pHist;
if isfield(par, 'pIn') && isfield(par, 'pOut')
    g1D.pIn = par.pIn;
    g1D.pOut = par.pOut;
else
    g1D.pIn = 0;
    g1D.pOut = 0;
end
end

function g1D = add_global_1d_blank_solid_pressure_view(out, par)
g1D = out.global1D;
if isstruct(g1D) && isfield(g1D, 'pHist')
    g1D.pBlanked = g1D.pHist;
    if isfield(out, 'deltaEHist') && isfield(out, 'deltaLHist')
        if isfield(par, 'gapFloor') && ~isempty(par.gapFloor)
            gap_floor_val = par.gapFloor;
        else
            gap_floor_val = 1.0e-8; % Default 10 nm contact threshold
        end
        contact_mask = (out.deltaEHist - out.deltaLHist) <= gap_floor_val;
        g1D.pBlanked(contact_mask) = NaN;
    end
end
end
% ========================================================================
% HELPER FUNCTIONS
% ========================================================================


function mesh = relax_surface_mesh_nodes(mesh, factor)
if nargin < 2, factor = 0.1; end
nodes = mesh.nodes;
conn = mesh.conn;
nNodes = size(nodes, 1);

adj = sparse(nNodes, nNodes);
for e = 1:size(conn, 1)
    nodes_e = conn(e, :);
    adj(nodes_e, nodes_e) = 1;
end
adj = adj - diag(diag(adj));

for i = 1:nNodes
    neighbors = find(adj(i, :));
    if ~isempty(neighbors)
        nodes(i, :) = (1 - factor) * nodes(i, :) + factor * mean(nodes(neighbors, :), 1);
    end
end
mesh.nodes = nodes;
end

function [R, Z, Uz] = velocity_field_for_plot(z, fluid, state, par)
if isfield(fluid, 'uzC') && ~isempty(fluid.uzC) && isfield(fluid, 'meshF')
    if isfield(fluid.meshF, 'Rp') && isfield(fluid.meshF, 'Zp')
        R = fluid.meshF.Rp;
        Z = fluid.meshF.Zp;
        Uz = fluid.uzC;
        return;
    elseif isfield(fluid.meshF, 'Ruz') && isfield(fluid.meshF, 'Zuz') && isfield(fluid.meshF, 'uz2D')
        R = fluid.meshF.Ruz;
        Z = fluid.meshF.Zuz;
        Uz = fluid.uz2D;
        return;
    end
end

[R, Z, Uz] = build_velocity_field(z, fluid.p, state.deltaL, ...
    state.deltaE, state.UwL, state.UwE, par);
end

function print_step_diagnostics(diag)
fprintf(['   diagnostics: fluidRes=%.3e, globalMass=%.3e, ', ...
    'dV/V=%.3e, fluxJump=%.3e, sourceInt=%.3e\n'], ...
    diag.fluidResidualInf, diag.globalMassResidual, ...
    diag.volumeChangeRel, diag.fluxJump, diag.sourceIntegral);
if isfield(diag, 'minJE')
    fprintf('   geometry: JEmin=%.3e, rEmin=%.3e m', ...
        diag.minJE, diag.minRadiusE);
    if isfield(diag, 'minJL') && isfinite(diag.minJL)
        fprintf(', JLmin=%.3e, rLmin=%.3e m', ...
            diag.minJL, diag.minRadiusL);
    end
    fprintf('\n');
end
end

function tf = should_retry_time_step(reason, dtAttempt, retryCount, par)
tf = false;
if ~isfield(par, 'enableAdaptiveTimeStep') || ~par.enableAdaptiveTimeStep
    return;
end

maxRetries = 6;
if isfield(par, 'maxTimeStepRetries') && isfinite(par.maxTimeStepRetries)
    maxRetries = par.maxTimeStepRetries;
end
if retryCount >= maxRetries
    return;
end

dtMin = 0;
if isfield(par, 'dtMin') && isfinite(par.dtMin)
    dtMin = par.dtMin;
end
if dtAttempt <= dtMin * (1 + 10*eps)
    return;
end

retryTokens = { ...
    'Negative or zero J', ...
    'Element inverted', ...
    'Non-positive radius', ...
    'Solid geometry guard failed', ...
    'Gap violates minGap', ...
    'violates minGap', ...
    'Fluid solve failed', ...
    'Pressure jump too large', ...
    'fsolve failed', ...
    'did not reach equilibrium', ...
    'residual too large', ...
    'line search failed', ...
    'did not converge'};

for k = 1:numel(retryTokens)
    if contains(reason, retryTokens{k})
        tf = true;
        return;
    end
end
end

function summary = summarize_time_step_adaptation(dtHist, retryHist, par, stoppedEarly, stopReason)
valid = isfinite(dtHist) & dtHist > 0;
summary = struct();
summary.enabled = isfield(par, 'enableAdaptiveTimeStep') && par.enableAdaptiveTimeStep;
summary.acceptedSteps = nnz(valid);
summary.totalRetries = sum(retryHist(valid));
summary.retrySteps = nnz(retryHist(valid) > 0);
summary.stoppedEarly = stoppedEarly;
summary.stopReason = stopReason;

if any(valid)
    summary.minDt = min(dtHist(valid));
    summary.maxDt = max(dtHist(valid));
    summary.finalDt = dtHist(find(valid, 1, 'last'));
    summary.maxRetriesInStep = max(retryHist(valid));
else
    summary.minDt = NaN;
    summary.maxDt = NaN;
    summary.finalDt = NaN;
    summary.maxRetriesInStep = 0;
end
end

function print_adaptive_summary(summary)
fprintf(['Adaptive stepping summary: accepted=%d, retrySteps=%d, ', ...
    'totalRetries=%d, minDt=%.3e, finalDt=%.3e\n'], ...
    summary.acceptedSteps, summary.retrySteps, summary.totalRetries, ...
    summary.minDt, summary.finalDt);
if summary.stoppedEarly
    fprintf('   stopped early: %s\n', compact_failure_reason(summary.stopReason));
end
end

function s = compact_failure_reason(reason)
s = regexprep(char(reason), '\s+', ' ');
maxChars = 180;
if numel(s) > maxChars
    s = [s(1:maxChars), '...'];
end
end

function state = attach_state_geometry_checks(state, meshE, meshL, par)
qE = solid_geometry_quality(meshE, state.uE, 'endothelium');
assert_solid_geometry_ok(qE, par);
state.geometryE = qE;

if ~isempty(meshL) && isfield(state, 'uL') && ~isempty(state.uL)
    qL = solid_geometry_quality(meshL, state.uL, 'leukocyte');
    assert_solid_geometry_ok(qL, par);
    state.geometryL = qL;
end
end

function warn_step_diagnostics(diag, par)
if ~isfield(par, 'diagnosticsEnabled') || ~par.diagnosticsEnabled
    return;
end

fluidTol = inf;
if isfield(par, 'diagnosticsWarnFluidResidual')
    fluidTol = par.diagnosticsWarnFluidResidual;
end
if isfinite(fluidTol) && diag.fluidResidualInf > fluidTol
    warning('Fluid diagnostic residual %.3e exceeds %.3e.', ...
        diag.fluidResidualInf, fluidTol);
end

volumeTol = inf;
if isfield(par, 'diagnosticsWarnVolumeJumpRel')
    volumeTol = par.diagnosticsWarnVolumeJumpRel;
end
if isfinite(volumeTol) && diag.volumeChangeRel > volumeTol
    warning('Relative fluid-volume change %.3e exceeds %.3e.', ...
        diag.volumeChangeRel, volumeTol);
end
end

function parL = leukocyte_solid_parameters(par)
parL = par;

if isfield(par, 'EL')
    parL.Ee = par.EL;
end
if isfield(par, 'nuL')
    parL.nuE = par.nuL;
end
if isfield(par, 'GL')
    parL.Ge = par.GL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ge = parL.Ee/(2*(1+parL.nuE));
end
if isfield(par, 'KL')
    parL.Ke = par.KL;
elseif isfield(parL, 'Ee') && isfield(parL, 'nuE')
    parL.Ke = parL.Ee/(3*(1-2*parL.nuE));
end
if isfield(par, 'etaL')
    parL.etaE = par.etaL;
end
if isfield(par, 'etaBulkL')
    parL.etaBulkE = par.etaBulkL;
end
if isfield(par, 'useViscoelasticLeukocyte')
    parL.useViscoelasticEndothelium = par.useViscoelasticLeukocyte;
end
% =========================================================================
% FORCE ACTIVE TRANSLOCATION PARAMETER PROPAGATION TO parL
% =========================================================================
if isfield(par, 'useActiveTranslocation')
    parL.useActiveTranslocation = par.useActiveTranslocation;
end
if isfield(par, 'fz_active_translocation')
    parL.fz_active_translocation = par.fz_active_translocation;
end
if isfield(par, 'z_pore_bottom')
    parL.z_pore_bottom = par.z_pore_bottom;
end
if isfield(par, 'sigma_head_decay')
    parL.sigma_head_decay = par.sigma_head_decay;
end
% =========================================================================
end

function warn_leukocyte_prestress_load_mismatch(SL, par)
if isfield(par, 'warnPrestressLoadMismatch') && ~par.warnPrestressLoadMismatch
    return;
end

[isMismatch, prestressLoad, runtimeLoad] = leukocyte_prestress_load_mismatch(SL, par);
if ~isMismatch
    return;
end

warning(['Leukocyte prestress load scale is %.3g Pa, while runtime ', ...
    'pIn/pOut are %.3g/%.3g Pa. Make sure the coupled initial ', ...
    'fluid/solid load is intentional, otherwise the first step may ', ...
    'mostly relax a prestress mismatch.'], ...
    prestressLoad, par.pIn, par.pOut);
end

function [isMismatch, prestressLoad, runtimeLoad] = leukocyte_prestress_load_mismatch(SL, par)
prestressLoad = 0;
if isfield(SL, 'P0') && isnumeric(SL.P0)
    p0Vals = SL.P0(:);
    p0Vals = p0Vals(isfinite(p0Vals));
    if ~isempty(p0Vals)
        prestressLoad = max(prestressLoad, max(abs(p0Vals)));
    end
end
if isfield(SL, 'trL') && isfield(SL.trL, 'normal') && isnumeric(SL.trL.normal)
    normalVals = SL.trL.normal(:);
    normalVals = normalVals(isfinite(normalVals));
    if ~isempty(normalVals)
        prestressLoad = max(prestressLoad, max(abs(normalVals)));
    end
end

runtimeLoad = max(abs([par.pIn, par.pOut]));
isMismatch = prestressLoad > max(10 * runtimeLoad, 1e-9);
end

function state = initial_state(z, par, meshE, uE_pre, deltaE_pre, ...
    meshL, interfaceL, uL_pre, deltaL_pre)
state = struct();

zVec = z(:);

state.UwL = zeros(size(zVec));
state.UwE = zeros(size(zVec));

state.zE = zVec;

if isfield(par, 'noLeukocyte') && par.noLeukocyte
    state.deltaL = zeros(size(zVec));
    state.zL = [];
    state.uL = [];
    state.uLPrev = [];
else
    state.uL = uL_pre(:);
    state.uLPrev = uL_pre(:);

    if use_RLout_fluid_interface_for_solid_leukocyte(par)
        state.deltaL = par.RLout * ones(size(zVec));
        state.zL = zVec;
        state.UwL = zeros(size(zVec));
    else
        state.deltaL = deltaL_pre(:);
        state.zL = zVec;
    end
end

state.deltaE = deltaE_pre(:);
state.p = linspace(par.pIn, par.pOut, numel(zVec)).';
state.pReduced = state.p;
state.p2D = state.p;

state.uE = uE_pre;
state.uEPrev = uE_pre;
state.pPrev = state.p;
state.dtPrev = par.dt;
end

function val = if_isfield_default(s, fieldname, defaultVal)
if isstruct(s) && isfield(s, fieldname)
    val = s.(fieldname);
else
    val = defaultVal;
end
end




