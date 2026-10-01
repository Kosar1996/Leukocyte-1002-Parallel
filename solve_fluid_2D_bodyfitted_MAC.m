function [fluid, ok, stopReason, meshF] = solve_fluid_2D_bodyfitted_MAC(z, old, state, par)
%SOLVE_FLUID_2D_BODYFITTED_MAC
% Body-fitted MAC wrapper with exact quadratic boundary flux formulation,
% mid-step midpoint mesh evaluation, and comprehensive diagnostic reporting.

ok = true;
stopReason = '';
meshF = [];

zOut = z(:);
NzOut = numel(zOut);

if NzOut < 2
    ok = false;
    stopReason = 'Body-fitted MAC Stokes fluid requires at least two axial points.';
    fluid = empty_fluid_2D_return(NzOut, nan(NzOut,1));
    return;
end

rlOut = state.deltaL(:);
reOut = state.deltaE(:);

useRLoutInner = use_RLout_fluid_interface_for_solid_leukocyte(par);
if useRLoutInner
    rlOut = par.RLout * ones(NzOut,1);
       fprintf('  if Detected Leukocyte r range softlube: [%.6e, %.6e] m\n', rlOut(3), rlOut(4));

end


% --- DYNAMICALLY DETERMINE ENDOTHELIUM AXIAL SPAN ---
if isfield(par, 'REout') && ~isempty(par.REout)
    baselineRE = par.REout;
else
    baselineRE = max(reOut);
end

endoNodes = abs(reOut - baselineRE) > 1.0e-7;
if any(endoNodes)
    zIndicesEndo = find(endoNodes);
    zMinEndo = zOut(zIndicesEndo(1));
    zMaxEndo = zOut(zIndicesEndo(end));
else
    zMinEndo = zOut(1);
    zMaxEndo = zOut(end);
end

outsideEndoSpan = (zOut < zMinEndo) | (zOut > zMaxEndo);
rlOut(outsideEndoSpan) = 0.0;

fprintf('\n=== [DYNAMIC ENDOTHELIUM SPAN AUDIT] ===\n');
fprintf('  Detected Endothelium Span : [%.6e, %.6e] m\n', zMinEndo, zMaxEndo);
fprintf('  Nodes Zeroed Outside Span : %d\n', sum(outsideEndoSpan));
fprintf('==========================================\n\n');

% --- FIX 1: Compute hOut prior to error guards ---
hOut = reOut - rlOut;

if numel(rlOut) ~= NzOut || numel(reOut) ~= NzOut
    ok = false;
    stopReason = 'Body-fitted MAC Stokes fluid: state.deltaL/deltaE size does not match z.';
    fluid = empty_fluid_2D_return(NzOut, hOut);
    return;
end

% --- RESTRICT MIN GAP CHECK STRICTLY WITHIN DYNAMIC ENDOTHELIUM SPAN ---
hOut_active = hOut((zOut >= zMinEndo) & (zOut <= zMaxEndo));
zOut_active = zOut((zOut >= zMinEndo) & (zOut <= zMaxEndo));

if any(hOut_active <= par.minGap)
    [hmin, imin_local] = min(hOut_active);
    imin = find(zOut == zOut_active(imin_local), 1);
    
    fprintf('\n================ [MIN GAP FAILURE AUDIT] ================\n');
    fprintf('  Failure at z = %.6e m, h = %.6e m (minGap = %.6e m)\n', zOut(imin), hmin, par.minGap);
    fprintf('  Endothelium radius reOut(z) = %.6e m\n', reOut(imin));
    fprintf('  Leukocyte radius   rlOut(z) = %.6e m\n', rlOut(imin));
    fprintf('  Dynamic Endothelium Span  : z in [%.6e, %.6e]\n', zMinEndo, zMaxEndo);
    fprintf('=========================================================\n\n');
    
    ok = false;
    stopReason = sprintf('Body-fitted MAC Stokes fluid: gap below minGap at z=%.6e, h=%.6e.', ...
        zOut(imin), hmin);
    fluid = empty_fluid_2D_return(NzOut, hOut);
    return;
end

zMin = par.zMin;
zMax = par.zMax;
if ~(isfinite(zMin) && isfinite(zMax) && zMax > zMin)
    zMin = zOut(1);
    zMax = zOut(end);
end

NzMAC = NzOut;
zF = linspace(zMin, zMax, NzMAC+1).';
zc = 0.5*(zF(1:end-1) + zF(2:end));

parMAC = par;
if isfield(par, 'NrFluid2D') && ~isempty(par.NrFluid2D)
    parMAC.Nr = par.NrFluid2D;
elseif ~isfield(parMAC, 'Nr') || isempty(parMAC.Nr)
    parMAC.Nr = 16;
end
parMAC.Nz = NzMAC;
parMAC.NrFluid2D = parMAC.Nr;

parMAC.pressurePenalty = 0;
if isfield(par, 'fluidPressurePenalty') && isfinite(par.fluidPressurePenalty)
    parMAC.pressurePenalty = par.fluidPressurePenalty;
end
if ~isfield(parMAC, 'pIn') || isempty(parMAC.pIn),  parMAC.pIn = 0;  end
if ~isfield(parMAC, 'pOut') || isempty(parMAC.pOut), parMAC.pOut = 0; end
parMAC.pInFun  = @(r,t) 0*r + parMAC.pIn;
parMAC.pOutFun = @(r,t) 0*r + parMAC.pOut;

oldDeltaL = old.deltaL(:);
oldDeltaE = old.deltaE(:);
if useRLoutInner
    oldDeltaL = par.RLout * ones(NzOut,1);
end

stateMAC = struct();
stateMAC.zc = zc;
stateMAC.zF = zF;
stateMAC.deltaL = safe_interp1_same_or_resample(zOut, rlOut, zc, 'state.deltaL');
stateMAC.deltaE = safe_interp1_same_or_resample(zOut, reOut, zc, 'state.deltaE');

oldMAC = struct();
oldMAC.deltaL = safe_interp1_same_or_resample(zOut, oldDeltaL, zc, 'old.deltaL');
oldMAC.deltaE = safe_interp1_same_or_resample(zOut, oldDeltaE, zc, 'old.deltaE');

if any(stateMAC.deltaE - stateMAC.deltaL <= par.minGap)
    hMAC = stateMAC.deltaE - stateMAC.deltaL;
    [hmin, imin] = min(hMAC);
    ok = false;
    stopReason = sprintf('Body-fitted MAC Stokes fluid: interpolated gap below minGap at z=%.6e, h=%.6e.', ...
        zc(imin), hmin);
    fluid = empty_fluid_2D_return(NzOut, hOut);
    return;
end

if useRLoutInner
    stateMAC.UwL = zeros(NzMAC,1);
elseif isfield(state, 'UwL') && numel(state.UwL) == NzOut
    stateMAC.UwL = safe_interp1_same_or_resample(zOut, state.UwL(:), zc, 'state.UwL');
else
    stateMAC.UwL = zeros(NzMAC,1);
end

if isfield(state, 'UwE') && numel(state.UwE) == NzOut
    stateMAC.UwE = safe_interp1_same_or_resample(zOut, state.UwE(:), zc, 'state.UwE');
else
    stateMAC.UwE = zeros(NzMAC,1);
end

bc = struct();
bc.uzL = stateMAC.UwL;
bc.uzE = stateMAC.UwE;

% =========================================================================
% DIAGNOSTIC AUDIT SUITE: FLUID VELOCITY CHECKS
% =========================================================================
fprintf('\n================== [FLUID DIAGNOSTIC AUDIT] ==================\n');

max_dDeltaE = max(abs(stateMAC.deltaE - oldMAC.deltaE));
max_dDeltaL = max(abs(stateMAC.deltaL - oldMAC.deltaL));
fprintf('[TEST 1: Interface Radius Update]\n');
fprintf('  max|deltaE^{n+1} - deltaE^n| = %.6e m\n', max_dDeltaE);
fprintf('  max|deltaL^{n+1} - deltaL^n| = %.6e m\n', max_dDeltaL);
if max_dDeltaE == 0 && max_dDeltaL == 0
    fprintf('  --> WARNING: Zero interface displacement passed from solid solver!\n');
else
    fprintf('  --> SUCCESS: Non-zero solid deformation detected.\n');
end

has_v_wall_E = isfield(par, 'v_wall_E') && numel(par.v_wall_E) == NzOut;
has_v_wall_L = isfield(par, 'v_wall_L') && numel(par.v_wall_L) == NzOut;
fprintf('[TEST 2: External Wall Velocity Override]\n');
fprintf('  par.v_wall_E field present: %s\n', string_on_off(has_v_wall_E));
fprintf('  par.v_wall_L field present: %s\n', string_on_off(has_v_wall_L));

deltaE_mid = 0.5 * (stateMAC.deltaE + oldMAC.deltaE);
deltaL_mid = 0.5 * (stateMAC.deltaL + oldMAC.deltaL);

urE_quadratic = (stateMAC.deltaE.^2 - oldMAC.deltaE.^2) ./ (2 * deltaE_mid * par.dt);
urE_linear    = (stateMAC.deltaE - oldMAC.deltaE) ./ par.dt;
urL_quadratic = (stateMAC.deltaL.^2 - oldMAC.deltaL.^2) ./ (2 * deltaL_mid * par.dt);

if has_v_wall_E
    graphUrE = safe_interp1_same_or_resample(zOut, par.v_wall_E(:), zc, 'par.v_wall_E');
    fprintf('  --> OVERRIDE ACTIVE: Endothelium using external par.v_wall_E.\n');
else
    graphUrE = urE_quadratic;
    fprintf('  --> QUADRATIC ACTIVE: Endothelium using exact volume flux (deltaE^2 - old^2)/(2*deltaE_mid*dt).\n');
end

if has_v_wall_L
    graphUrL = safe_interp1_same_or_resample(zOut, par.v_wall_L(:), zc, 'par.v_wall_L');
    fprintf('  --> OVERRIDE ACTIVE: Leukocyte using external par.v_wall_L.\n');
else
    graphUrL = urL_quadratic;
    fprintf('  --> QUADRATIC ACTIVE: Leukocyte using exact volume flux (deltaL^2 - old^2)/(2*deltaL_mid*dt).\n');
end

fprintf('[TEST 3: Quadratic Velocity Path]\n');
fprintf('  Endothelium max|ur_quadratic| = %.6e m/s\n', max(abs(urE_quadratic)));
fprintf('  Endothelium max|ur_linear|    = %.6e m/s\n', max(abs(urE_linear)));
fprintf('  Quadratic vs Linear relative diff = %.2f%%\n', ...
    100 * max(abs(urE_quadratic - urE_linear)) / max(abs(urE_linear) + realmin));

useSlopeKinematics = isfield(par, 'useSlopeAwareWallKinematics') && ...
    par.useSlopeAwareWallKinematics;
fprintf('[TEST 4: Slope Kinematics Status]\n');
fprintf('  useSlopeAwareWallKinematics = %s\n', string_on_off(useSlopeKinematics));

if useSlopeKinematics
    slopeL = curve_slope_1d(stateMAC.zc, stateMAC.deltaL);
    slopeE = curve_slope_1d(stateMAC.zc, stateMAC.deltaE);

    maxSlope = 0.75;
    slopeL = sign(slopeL) .* min(abs(slopeL), maxSlope);
    slopeE = sign(slopeE) .* min(abs(slopeE), maxSlope);

    uzL_rel = bc.uzL;
    uzE_rel = bc.uzE;

    if isfield(oldMAC, 'UwL') && numel(oldMAC.UwL) == NzMAC
        uzL_mesh = 0.5 * (bc.uzL + oldMAC.UwL);
        uzL_rel  = bc.uzL - uzL_mesh;
    end
    if isfield(oldMAC, 'UwE') && numel(oldMAC.UwE) == NzMAC
        uzE_mesh = 0.5 * (bc.uzE + oldMAC.UwE);
        uzE_rel  = bc.uzE - uzE_mesh;
    end

    bc.urL = graphUrL + slopeL .* uzL_rel;
    bc.urE = graphUrE + slopeE .* uzE_rel;
    fprintf('  Max surface slope (dr/dz) = %.4e | Max relative uz = %.4e m/s\n', ...
        max(abs(slopeE)), max(abs(uzE_rel)));
else
    bc.urL = graphUrL;
    bc.urE = graphUrE;
end
fprintf('==============================================================\n\n');

try
    stateMAC_mid = stateMAC;
    stateMAC_mid.deltaE = deltaE_mid;
    stateMAC_mid.deltaL = deltaL_mid;

    meshMAC = build_body_fitted_gap_mesh(stateMAC_mid, parMAC);
    tNow = 0;
    if isfield(state, 't') && isfinite(state.t)
        tNow = state.t;
    end

    if isfield(par, 'useUnsteadyStokes') && par.useUnsteadyStokes
        if isfield(old, 'ur2DFaceField') && ~isempty(old.ur2DFaceField) && ...
                isequal(size(old.ur2DFaceField), [parMAC.Nr+1, parMAC.Nz]) && ...
                isfield(old, 'Rur2DField') && ~isempty(old.Rur2DField) && ...
                isequal(size(old.Rur2DField), [parMAC.Nr+1, parMAC.Nz])
            bc.urPrev = interp_facefield_radial_to_mesh(old.ur2DFaceField, old.Rur2DField, meshMAC.Rur);
        else
            bc.urPrev = zeros(parMAC.Nr+1, parMAC.Nz);
        end
        if isfield(old, 'uz2DFaceField') && ~isempty(old.uz2DFaceField) && ...
                isequal(size(old.uz2DFaceField), [parMAC.Nr, parMAC.Nz+1]) && ...
                isfield(old, 'Ruz2DField') && ~isempty(old.Ruz2DField) && ...
                isequal(size(old.Ruz2DField), [parMAC.Nr, parMAC.Nz+1])
            bc.uzPrev = interp_facefield_radial_to_mesh(old.uz2DFaceField, old.Ruz2DField, meshMAC.Ruz);
        else
            bc.uzPrev = zeros(parMAC.Nr, parMAC.Nz+1);
        end
    end

    fluidMAC = solve_stokes_bodyfitted_MAC(meshMAC, stateMAC, bc, parMAC, tNow);
    fluidMAC.bc = bc;
    [tauLMAC, tauEMAC] = estimate_wall_shear_bodyfitted(meshMAC, fluidMAC, stateMAC, parMAC);
catch ME
    ok = false;
    stopReason = ['Body-fitted MAC Stokes fluid failed: ', ME.message];
    fluid = empty_fluid_2D_return(NzOut, hOut);
    return;
end

pMidMAC = extract_midgap_pressure_line_bodyfitted(fluidMAC);
pLMAC = fluidMAC.P(1,:).';
pEMAC = fluidMAC.P(end,:).';
uzLMAC = fluidMAC.uz(1,:).';
uzEMAC = fluidMAC.uz(end,:).';

pOutVec    = safe_interp1_same_or_resample(zc, pMidMAC(:), zOut, 'pMidMAC');
pLOutVec   = safe_interp1_same_or_resample(zc, pLMAC(:),    zOut, 'pLMAC');
pEOutVec   = safe_interp1_same_or_resample(zc, pEMAC(:),    zOut, 'pEMAC');

pOutVec(1)   = parMAC.pIn;
pOutVec(end) = parMAC.pOut;
pLOutVec(1)   = parMAC.pIn;
pLOutVec(end) = parMAC.pOut;
pEOutVec(1)   = parMAC.pIn;
pEOutVec(end) = parMAC.pOut;

tauLOutVec = safe_interp1_same_or_resample(zc, tauLMAC(:),  zOut, 'tauLMAC');
tauEOutVec = safe_interp1_same_or_resample(zc, tauEMAC(:),  zOut, 'tauEMAC');
uzLOutVec  = safe_interp1_same_or_resample(zc, uzLMAC(:),   zOut, 'uzLOutVec');
uzEOutVec  = safe_interp1_same_or_resample(zc, uzEMAC(:),   zOut, 'uzEOutVec');

Qfaces = compute_bodyfitted_axisym_flux_faces(meshMAC, fluidMAC);
zQfaces = meshMAC.zF(:);
zQout = 0.5*(zOut(1:end-1) + zOut(2:end));
Qout = safe_interp1_same_or_resample(zQfaces, Qfaces(:), zQout, 'Qfaces');

fluid = struct();
fluid.meshType = 'bodyfitted_MAC';
fluid.leukocyteInnerBoundary = ternary_local(useRLoutInner, 'RLout_reference', 'state_deltaL');
fluid.usesRLoutFluidInterfaceForSolidLeukocyte = useRLoutInner;
fluid.meshF = meshMAC;
meshF = meshMAC;

fluid.p = pOutVec(:);
fluid.pL = pLOutVec(:);
fluid.pE = pEOutVec(:);
fluid.Q = Qout(:);
fluid.tauL = tauLOutVec(:);
fluid.tauE = tauEOutVec(:);
fluid.uzL = uzLOutVec(:);
fluid.uzE = uzEOutVec(:);
fluid.gap = hOut(:);

fluid.P = fluidMAC.P;
fluid.ur = fluidMAC.ur;
fluid.uz = fluidMAC.uz;
fluid.urC = fluidMAC.urC;
fluid.uzC = fluidMAC.uzC;
fluid.divInf = fluidMAC.divInf;
fluid.linRes = fluidMAC.linRes;
fluid.condNumber = fluidMAC.condNumber;
fluid.Nunknown = fluidMAC.Nunknown;
fluid.bc = bc;
[fluid.tractionL, fluid.tractionE] = compute_bodyfitted_wall_traction(meshMAC, fluidMAC, parMAC);
fluid.pAvg = radial_average_pressure_bodyfitted(fluidMAC.P, meshMAC);
fluid.pMid = pMidMAC(:);

fluid.ur2D = fluidMAC.urC(:);
fluid.uz2D = fluidMAC.uzC(:);
fluid.pCell = fluidMAC.P(:);
fluid.sigmaCell = [];
fluid.cellCenter = [meshMAC.Rp(:), meshMAC.Zp(:)];
end

function str = string_on_off(val)
if islogical(val)
    if val, str = 'ON'; else, str = 'OFF'; end
elseif isnumeric(val)
    if val ~= 0, str = 'ON'; else, str = 'OFF'; end
else
    str = char(val);
end
end