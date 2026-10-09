function stateNew = solve_monolithic_two_solids_fsolve_timestep( ...
    old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par)
%SOLVE_MONOLITHIC_TWO_SOLIDS_FSOLVE_TIMESTEP
% Residual-only fully monolithic coupling for deformable leukocyte:
% unknown vector y = [uE_free/uScale; uL_free/uScale; p_internal/pScale].

% Ensure 1-based indexing for boundary and support arrays
if any(baseE(:) == 0), baseE = baseE + 1; end
if any(baseL(:) == 0), baseL = baseL + 1; end
if exist('interfaceE', 'var') && any(interfaceE(:) == 0), interfaceE = interfaceE + 1; end

ndofE = size(meshE.nodes,1) * 2;
ndofL = size(meshL.nodes,1) * 2;
N = numel(z);

[fixE, valsE] = solid_support_conditions(baseE, par.supportE);
fixE = unique(fixE(:));
valsE = valsE(:);
freeE = setdiff((1:ndofE).', fixE);

[fixL, valsL] = solid_support_conditions(baseL, par.supportL);
fixL = unique(fixL(:));
valsL = valsL(:);
freeL = setdiff((1:ndofL).', fixL);

if isfield(par, 'twoSolidInterfaceOnlyTest') && par.twoSolidInterfaceOnlyTest
    interfaceDofsE = reshape([2*interfaceE(:)-1, 2*interfaceE(:)].', [], 1);
    interfaceDofsL = reshape([2*interfaceL(:)-1, 2*interfaceL(:)].', [], 1);
    freeE = intersect(freeE, interfaceDofsE);
    freeL = intersect(freeL, interfaceDofsL);
    fprintf('FAST TEST: using interface-only solid DOFs: nE=%d, nL=%d, p=%d\n', ...
        numel(freeE), numel(freeL), max(numel(z)-2,0));
end

% Synchronized monolithic scaling factors
JuE = par.uScaleMono;
JuL = par.uScaleMono;
Jp  = par.pScaleMono;

[solidTarget, fluidTarget] = monolithic_residual_targets(par);
solidTargetE = solidTarget;
solidTargetL = solidTarget;

if isfield(par, 'monoSolidAbsTol') && isfinite(par.monoSolidAbsTol) && par.monoSolidAbsTol > 0
    solidTargetE = par.monoSolidAbsTol;
    solidTargetL = par.monoSolidAbsTol;
end
if isfield(par, 'monoFluidAbsTol') && isfinite(par.monoFluidAbsTol) && par.monoFluidAbsTol > 0
    fluidTarget = par.monoFluidAbsTol;
end

% =========================================================================
% MANDATORY ACTIVE TRANSLOCATION PARAMETER PROPAGATION TO parL only when parL doesn't have its value
% =========================================================================
if isstruct(par) && isfield(par, 'useActiveTranslocation') && par.useActiveTranslocation
    parL.useActiveTranslocation  = true;
    if ~isfield(parL, 'fz_active_translocation')&&isfield(par, 'fz_active_translocation')
        parL.fz_active_translocation = par.fz_active_translocation;
    end
    if isfield(par, 'z_pore_bottom')
        parL.z_pore_bottom = par.z_pore_bottom;
    end
    if isfield(par, 'sigma_head_decay')
        parL.sigma_head_decay = par.sigma_head_decay;
    end
end

% Boundary condition enforcement for base state
uE0 = old.uE;
uL0 = old.uL;
p0  = old.p;
if ~all(isfinite(p0)) || numel(p0) ~= N
    p0 = linspace(par.pIn, par.pOut, N).';
end
uE0(fixE) = valsE;
uL0(fixL) = valsL;
p0(1)   = par.pIn;
p0(end) = par.pOut;

% =========================================================================
% MONOLITHIC PREDICTOR WARM-START INITIAL GUESS VECTOR (y0)
% =========================================================================
hasEPrev = isfield(old, 'uEPrev') && ~isempty(old.uEPrev) && ~isequal(old.uE, old.uEPrev);
hasLPrev = isfield(old, 'uLPrev') && ~isempty(old.uLPrev) && ~isequal(old.uL, old.uLPrev);

usePredictor = isfield(par, 'useMonoPredictor') && par.useMonoPredictor && hasEPrev && hasLPrev;

if usePredictor
    predScale = 1.0;
    if isfield(old, 'dtPrev') && old.dtPrev > 0
        predScale = min(1.0, par.dt / old.dtPrev);
    end

    % Extrapolate displacements u^{14}_pred = u^{13} + 1.0 * (u^{13} - u^{12})
    uEPred = old.uE + predScale * (old.uE - old.uEPrev);
    uLPred = old.uL + predScale * (old.uL - old.uLPrev);

    pPred = p0;
    uEPred(fixE) = valsE;
    uLPred(fixL) = valsL;

    y0 = [
        uEPred(freeE) / JuE
        uLPred(freeL) / JuL
        pPred(2:end-1) / Jp
    ];

    fprintf('   [fsolve Predictor] Active: Extrapolated displacements y0 with base p0.\n');
    fprintf('                      ||uEPred-uE0|| = %.3e m | ||uLPred-uL0|| = %.3e m\n', ...
        norm(uEPred - uE0), norm(uLPred - uL0));
else
    y0 = [
        uE0(freeE) / JuE
        uL0(freeL) / JuL
        p0(2:end-1) / Jp
    ];
    fprintf('   [fsolve Predictor] Fallback: y0 initialized without displacement history.\n');
end

% =========================================================================
% FSOLVE ACTIVE FORCE TWO-PASS RESIDUAL AUDIT (EXECUTED AFTER y0 IS SET)
% =========================================================================
fprintf('\n=== [FSOLVE ACTIVE FORCE CONTRIBUTION AUDIT] ===\n');
fprintf('  - parL.useActiveTranslocation     : %d\n', ...
    isfield(parL, 'useActiveTranslocation') && parL.useActiveTranslocation);

if isfield(parL, 'fz_active_translocation')
    fprintf('  - parL.fz_active_translocation    : %.3e N/m^3\n', parL.fz_active_translocation);
else
    fprintf('  - parL.fz_active_translocation    : 0.000e+00 N/m^3\n');
end

[uE_diag, uL_diag, p_diag] = unpack_two_solid_y( ...
    y0, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);

% Pass 1: Compute residual WITH active force
[~, RL_active, ~] = monolithic_two_solids_residual_unscaled( ...
    uE_diag, uL_diag, p_diag, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

% Pass 2: Compute residual WITHOUT active force (passive baseline)
% Pass 2: Compute true passive residual baseline for diagnostic audit
par_passive = par;
par_passive.useActiveTranslocation = false;

parL_passive = parL;
parL_passive.useActiveTranslocation = false;
parL_passive.fz_active_translocation = 0.0;

[~, RL_passive, ~] = monolithic_two_solids_residual_unscaled( ...
    uE_diag, uL_diag, p_diag, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL_passive, z, par_passive, freeE, freeL);

f_active_extracted = RL_passive - RL_active;
norm_f_active      = norm(f_active_extracted);
max_f_active       = max(abs(f_active_extracted));
norm_RL_passive    = norm(RL_passive);

fprintf('  - ||RL_passive|| (Elastic+Viscous+Fluid) : %.6e N\n', norm_RL_passive);
fprintf('  - ||RL_active||  (Total Residual)         : %.6e N\n', norm(RL_active));
fprintf('  - ||f_active(freeL)||                     : %.6e N\n', norm_f_active);
fprintf('  - Max Single DOF Active Force             : %.6e N\n', max_f_active);

if norm_RL_passive > 0
    fprintf('  - Active Force Weight Ratio               : %.2f%%\n', ...
        (norm_f_active / norm_RL_passive) * 100);
end

if norm_f_active == 0
    fprintf('  --> WARNING: Active force vector is ZERO! Check z_pore_bottom and element coordinates.\n');
else
    fprintf('  --> SUCCESS: Active force IS non-zero and actively modifying fsolve residual.\n');
end
fprintf('================================================\n\n');

nE = numel(freeE);
nL = numel(freeL);

% =========================================================================
% HARD-OVERRIDE: GUARANTEE parL ACTIVE TRANSLOCATION BEFORE FSOLVE HANDLE
% =========================================================================
if (isstruct(par) && isfield(par, 'useActiveTranslocation') && par.useActiveTranslocation) || ...
   (isstruct(parL) && isfield(parL, 'useActiveTranslocation') && parL.useActiveTranslocation)
    
    parL.useActiveTranslocation  = true;
    if ~isfield(parL, 'fz_active_translocation')&&isfield(par, 'fz_active_translocation')
        parL.fz_active_translocation = par.fz_active_translocation;
    end
    if isfield(par, 'z_pore_bottom')
        parL.z_pore_bottom = par.z_pore_bottom;
    end
    if isfield(par, 'sigma_head_decay')
        parL.sigma_head_decay = par.sigma_head_decay;
    end
end

% % =========================================================================
% % HARD OVERRIDE: FORCE FINITE-DIFFERENCE FSOLVE JACOBIAN FOR ACTIVE FORCE
% % =========================================================================
% par.useTwoSolidSemiAnalyticalJacobian = false;
% % =========================================================================

useSemiJac = isfield(par, 'useTwoSolidSemiAnalyticalJacobian') && ...
    par.useTwoSolidSemiAnalyticalJacobian;


if useSemiJac
    fprintf('jacobian scaled\n')
    fun = @(y) monolithic_two_solids_residual_jacobian_scaled( ...
        y, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, ...
        z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
        JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget);
else
fprintf('scaled\n')
    fun = @(y) monolithic_two_solids_residual_scaled( ...
        y, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, ...
        z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
        JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget);
end

if useSemiJac
    [R_y0, ~] = fun(y0);
else
    R_y0 = fun(y0);
end
fprintf('   [y0 Residual Verification] Initial norm ||f(y0)||^2 = %.5e\n\n', norm(R_y0)^2);

fsolveDisplay = 'iter';
if isfield(par, 'fsolveDisplay')
    fsolveDisplay = par.fsolveDisplay;
end

maxFun = max(800, 5 * max(numel(y0), 1));
if isfield(par, 'maxFunctionEvaluationsTwoSolid') && isfinite(par.maxFunctionEvaluationsTwoSolid)
    maxFun = par.maxFunctionEvaluationsTwoSolid;
end

if isfield(par, 'checkTwoSolidAnalyticalJacobian') && par.checkTwoSolidAnalyticalJacobian
    check_two_solid_jacobian_columns(fun, y0, nE, nL, N, par);
    par.checkTwoSolidAnalyticalJacobian = false;
end

bestY = y0;
bestR = [];
bestExitflag = NaN;
bestOutput = struct('iterations', 0);
bestScaledRes = inf;

algorithms = {'trust-region-dogleg', 'levenberg-marquardt'};

for attempt = 1:numel(algorithms)
    typicalX_scale = ones(size(y0));

    opts = optimoptions('fsolve', ...
        'Algorithm', algorithms{attempt}, ...
        'Display', fsolveDisplay, ...
        'SpecifyObjectiveGradient', useSemiJac, ...
        'ScaleProblem', 'Jacobian', ...
        'TypicalX', typicalX_scale, ...
        'FunctionTolerance', 1e-6, ...
        'StepTolerance', 1e-6, ...
        'OptimalityTolerance', 1e-6, ...
        'MaxIterations', par.maxNewtonMono, ...
        'MaxFunctionEvaluations', maxFun);

    [yAttempt, RAttempt, exitflagAttempt, outputAttempt] = fsolve(fun, y0, opts);
    scaledAttempt = norm(RAttempt, inf);

    if scaledAttempt < bestScaledRes || (exitflagAttempt > 0 && bestExitflag <= 0)
        bestY = yAttempt;
        bestR = RAttempt;
        bestExitflag = exitflagAttempt;
        bestOutput = outputAttempt;
        bestScaledRes = scaledAttempt;
    end
    if exitflagAttempt > 0
        break;
    end
end

ySol     = bestY;
Rsol     = bestR;
exitflag = bestExitflag;
output   = bestOutput;

[uESol, uLSol, ~] = unpack_two_solid_y( ...
    ySol, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
assert_solid_geometry_ok( ...
    solid_geometry_quality(meshE, uESol, 'endothelium fsolve solution'), par);
assert_solid_geometry_ok( ...
    solid_geometry_quality(meshL, uLSol, 'leukocyte fsolve solution'), par);

[solidENorm, solidLNorm, fluidNorm, gapMin] = two_solid_residual_norms_from_y( ...
    ySol, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);

fprintf('   [fsolve Mono Sub-Step] Exitflag = %d | Iters = %d | ScaledRes = %.3e | minGap = %.3e um\n', ...
    exitflag, output.iterations, norm(Rsol, inf), gapMin * 1e6);

if gapMin <= par.minGap
    fprintf('   [WARNING: GAP VIOLATION] minGap threshold reached (h_min = %.4e um <= %.4e um).\n', ...
        gapMin * 1e6, par.minGap * 1e6);
end

scaledRes = norm(Rsol, inf);

acceptByPhysicalResidual = ...
    solidENorm < solidTargetE && solidLNorm < solidTargetL && fluidNorm < fluidTarget;
acceptByScaledResidual = scaledRes < 1e-4;

if exitflag <= 0 && ~(acceptByPhysicalResidual || acceptByScaledResidual)
    fprintf('\n   [MONOLITHIC GUARD] fsolve exitflag=%d (scaledRes = %.3e). Executing Armijo Line-Search...\n', ...
        exitflag, scaledRes);

    normR0 = norm(R_y0);
    dy = ySol - y0;
    alpha_ls = 0.50;
    ls_success = false;

    for ls_iter = 1:5
        yCandidate = y0 + alpha_ls * dy;
        if useSemiJac
            [RCand, ~] = fun(yCandidate);
        else
            RCand = fun(yCandidate);
        end

        % Require significant residual reduction (norm < 0.8 * normR0) to accept line search
        if norm(RCand) < 0.80 * normR0
            fprintf('   [Line Search] Accepted step size alpha = %.4f (Norm: %.3e -> %.3e)\n', ...
                alpha_ls, normR0, norm(RCand));
            ySol = yCandidate;
            ls_success = true;
            break;
        end
        alpha_ls = alpha_ls * 0.5;
    end

    % CHANGE 1 ENFORCEMENT:
    % Throw an explicit error if line search fails OR if the post-line-search residual remains large (> 1.0e-1).
    % This triggers the try-catch block in softlube_run_case_global_coupled.m and forces a dt time-step reduction.
    finalScaledRes = norm(fun(ySol), inf);
    if ~ls_success || finalScaledRes > 1.0e-1
        error('Monolithic:fsolveFailed', ...
            'fsolve failed to converge (exitflag = %d, finalScaledRes = %.3e). Forcing time-step reduction.', ...
            exitflag, finalScaledRes);
    end

    [uESol, uLSol, ~] = unpack_two_solid_y( ...
        ySol, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
end

if gapMin <= par.minGap
    error('FSI:GapViolation', 'Two-solid monolithic solution violates minGap. gapMin = %.6e m', gapMin);
end

if exitflag <= 0 && acceptByScaledResidual
    fprintf('   accepting fsolve result: scaledRes %.3e < 1e-4 strict full-DOF-test cutoff.\n', scaledRes);
end

stateNew = build_two_solid_state_from_y( ...
    ySol, old, meshE, interfaceE, meshL, interfaceL, z, par, ...
    freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp);

if isfield(par, 'dt') && par.dt > 0
    v_max_cap = 0.020;
    if isfield(par, 'v_max_cap') && isfinite(par.v_max_cap) && par.v_max_cap > 0
        v_max_cap = par.v_max_cap;
    end

    v_E_raw = (stateNew.deltaE - old.deltaE) / par.dt;
    if max(abs(v_E_raw)) > v_max_cap
        fprintf('   -> [Kinematic Clamp] Endothelium v_raw (%.3e m/s) clamped to %.3e m/s.\n', ...
            max(abs(v_E_raw)), v_max_cap);
        v_E_clamped = sign(v_E_raw) .* min(abs(v_E_raw), v_max_cap);
        stateNew.deltaE = old.deltaE + v_E_clamped * par.dt;
        stateNew.v_wall_E = v_E_clamped;
    end

    if isfield(stateNew, 'deltaL') && ~isempty(stateNew.deltaL)
        v_L_raw = (stateNew.deltaL - old.deltaL) / par.dt;
        if max(abs(v_L_raw)) > v_max_cap
            v_L_clamped = sign(v_L_raw) .* min(abs(v_L_raw), v_max_cap);
            stateNew.deltaL = old.deltaL + v_L_clamped * par.dt;
            stateNew.v_wall_L = v_L_clamped;
        end
    end
end
end

function check_two_solid_jacobian_columns(fun, y0, nE, nL, N, par)
fprintf('\nSelected-column check for two-solid analytical Jacobian:\n');
[R0, J0] = fun(y0);
cand = [];
labels = {};

if nE >= 1
    cand(end+1) = 1;
    labels{end+1} = 'uE first';
    cand(end+1) = max(1, round(nE/2));
    labels{end+1} = 'uE middle';
end
if nL >= 1
    cand(end+1) = nE + 1;
    labels{end+1} = 'uL first';
    cand(end+1) = nE + max(1, round(nL/2));
    labels{end+1} = 'uL middle';
end
nP = max(N-2,0);
if nP >= 1
    cand(end+1) = nE + nL + max(1, round(nP/2));
    labels{end+1} = 'p middle';
end

cand = unique(cand, 'stable');
for k = 1:numel(cand)
    j = cand(k);
    h = 1e-6 * max(1, abs(y0(j)));
    yp = y0; ym = y0;
    yp(j) = yp(j) + h;
    ym(j) = ym(j) - h;

    Rp = fun(yp);
    Rm = fun(ym);
    fdCol = (Rp - Rm) / (2*h);
    anCol = J0(:,j);

    absErr = norm(anCol - fdCol, inf);
    relErr = absErr / max([norm(fdCol, inf), norm(anCol, inf), eps]);
    fprintf('   %-10s col %6d: relErr = %.3e, absErr = %.3e\n', ...
        labels{k}, j, relErr, absErr);
end
fprintf('   initial scaled residual norm = %.3e\n\n', norm(R0, inf));
end

function [solidENorm, solidLNorm, fluidNorm, gapMin] = two_solid_residual_norms_from_y( ...
    y, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp)

if isstruct(par) && isfield(par, 'useActiveTranslocation') && par.useActiveTranslocation
    parL.useActiveTranslocation  = true;
    parL.fz_active_translocation = par.fz_active_translocation;
    parL.z_pore_bottom           = par.z_pore_bottom;
    parL.sigma_head_decay        = par.sigma_head_decay;
end

[uE, uL, p] = unpack_two_solid_y( ...
    y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);

[deltaE, ~] = monolithic_interface_kinematics_value_only( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, ~] = monolithic_interface_kinematics_value_only( ...
    meshL, uL, old.uL, interfaceL, z, par);
gapMin = min(deltaE - deltaL);

try
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshE, uE, 'endothelium fsolve candidate'), par);
    assert_solid_geometry_ok( ...
        solid_geometry_quality(meshL, uL, 'leukocyte fsolve candidate'), par);
catch
    solidENorm = inf;
    solidLNorm = inf;
    fluidNorm = inf;
    return;
end

if gapMin <= par.minGap
    solidENorm = inf;
    solidLNorm = inf;
    fluidNorm = inf;
    return;
end

try
    [RE, RL, RF] = monolithic_two_solids_residual_unscaled( ...
        uE, uL, p, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

    solidENorm = norm(RE, inf);
    solidLNorm = norm(RL, inf);
    fluidNorm  = norm(RF, inf);
catch
    solidENorm = inf;
    solidLNorm = inf;
    fluidNorm  = inf;
end
end

function stateNew = build_two_solid_state_from_y( ...
    y, old, meshE, interfaceE, meshL, interfaceL, z, par, ...
    freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp)

[uE, uL, p] = unpack_two_solid_y( ...
    y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);

[deltaE, UwE] = monolithic_interface_kinematics_value_only( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, UwL] = monolithic_interface_kinematics_value_only( ...
    meshL, uL, old.uL, interfaceL, z, par);

[Q, ~, ~, tauL, tauE, uzL, uzE] = ...
    local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

stateNew = old;
stateNew.uE = uE;
stateNew.uL = uL;
stateNew.uEPrev = old.uE;
stateNew.uLPrev = old.uL;

stateNew.deltaE = deltaE;
stateNew.deltaL = deltaL;
stateNew.UwE = UwE;
stateNew.UwL = UwL;

stateNew.p = p;
stateNew.pReduced = p;
stateNew.pPrev = old.p;
stateNew.dtPrev = par.dt;

stateNew.Q = Q;
stateNew.tauE = tauE;
stateNew.tauL = tauL;
stateNew.uzE = uzE;
stateNew.uzL = uzL;

if use_global2d_pressure_traction(par)
    [pLoadE, pLoadL, ~, ~, pressure2D] = global2d_pressure_traction_loads( ...
        z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);
    stateNew.pEGlobal2D = pLoadE;
    stateNew.pLGlobal2D = pLoadL;
    stateNew.global2DPressureTraction = pressure2D;
end
end
