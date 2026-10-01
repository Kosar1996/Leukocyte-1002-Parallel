% =========================================================================
% HEADER SUMMARY OF CHANGES:
% 1. STRICT DT SYNCHRONIZATION: Prioritize par.dt passed from the current solver
%    retry step over state.dt or fallback defaults to ensure moving-wall velocities 
%    (v_wall_E, v_wall_L) strictly match the actual time-step increment.
% 2. Safe shape assertions on interface displacement vectors to prevent shape 
%    mismatch during MAC boundary velocity injection.
% 3. Parameter fallbacks for smooth hybrid blending settings.
% 4. UNIFIED TRACTION & PRESSURE CAPPING: Enforce max pressure/shear traction 
%    caps on outgoing fluid struct to prevent sub-micron pressure runaway feedback loops.
% =========================================================================

function [fluid, ok, stopReason, meshF] = solve_selected_poststep_fluid(z, old, state, par)
%SOLVE_SELECTED_POSTSTEP_FLUID Dispatch the post-monolithic fluid evaluation.
%   Calculates boundary interface motion, dispatches to the requested solver mode,
%   and enforces traction/pressure guards prior to solid sub-solver transfer.

zVec = z(:);

% Ensure smooth hybrid blending fields have safe defaults if not specified
if ~isfield(par, 'useSmoothHybridBlending')
    par.useSmoothHybridBlending = true;
end
if ~isfield(par, 'hybridTransitionBuffer') || isempty(par.hybridTransitionBuffer)
    par.hybridTransitionBuffer = 0.3e-6;
end

% Set safe defaults for hydrodynamic traction caps if omitted
if ~isfield(par, 'maxFluidPressureCap') || isempty(par.maxFluidPressureCap)
    par.maxFluidPressureCap = 3000.0; % Pa (Upper bound for sub-micron squeezing)
end
if ~isfield(par, 'maxFluidShearCap') || isempty(par.maxFluidShearCap)
    par.maxFluidShearCap = 500.0; % Pa (Upper bound for interface shear traction)
end

% =========================================================================
% INTERFACE MOTION CALCULATION: Moving-Wall Boundary Velocity Injection
% =========================================================================
% STRICT PRIORITY: Use current trial step size par.dt directly from driver
dt = 1.0e-5;
if isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0
    dt = par.dt;
elseif isfield(state, 'dt') && isfinite(state.dt) && state.dt > 0
    dt = state.dt;
end

% Compute radial interface displacement velocities via backward difference
if isfield(state, 'deltaE') && isfield(old, 'deltaE') && ~isempty(old.deltaE) && ...
        numel(state.deltaE) == numel(zVec) && numel(old.deltaE) == numel(zVec)
    par.v_wall_E = (state.deltaE(:) - old.deltaE(:)) / dt;
else
    par.v_wall_E = zeros(size(zVec));
end

if isfield(state, 'deltaL') && isfield(old, 'deltaL') && ~isempty(old.deltaL) && ...
        numel(state.deltaL) == numel(zVec) && numel(old.deltaL) == numel(zVec)
    par.v_wall_L = (state.deltaL(:) - old.deltaL(:)) / dt;
else
    par.v_wall_L = zeros(size(zVec));
end

% Compute axial interface velocities if available
if isfield(state, 'UwE') && ~isempty(state.UwE) && numel(state.UwE) == numel(zVec)
    par.v_wall_zE = state.UwE(:);
else
    par.v_wall_zE = zeros(size(zVec));
end

if isfield(state, 'UwL') && ~isempty(state.UwL) && numel(state.UwL) == numel(zVec)
    par.v_wall_zL = state.UwL(:);
else
    par.v_wall_zL = zeros(size(zVec));
end
% =========================================================================

% STRICT PRIORITY: Full 2D Body-Fitted MAC Mode evaluates first
if isfield(par, 'useFull2DFluid') && par.useFull2DFluid && ...
   (~isfield(par, 'useHybridGap1DExterior2DFluid') || ~par.useHybridGap1DExterior2DFluid)
   
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_2D_bodyfitted_MAC(zVec, old, state, par);
        
elseif use_hybrid_gap1d_exterior2d_fluid(par)
    [fluid, ok, stopReason, meshF] = ...
        solve_fluid_hybrid_gap1d_exterior2d_withnodes(zVec, old, state, par);
else
    [fluid, ok, stopReason] = solve_fluid_reynolds_slip(zVec, old, state, par);
    meshF = [];
end

% =========================================================================
% UNIFIED TRACTION & PRESSURE CAPPING (SAFEGUARD AGAINST PRESSURE RUNAWAY)
% =========================================================================
if ok && isstruct(fluid)
    pCap   = par.maxFluidPressureCap;
    tauCap = par.maxFluidShearCap;

    % Diagnostic tracking
    if isfield(fluid, 'p') && ~isempty(fluid.p) && max(abs(fluid.p)) > pCap
        fprintf('    [Fluid Cap Active] Peak pressure (%.1f Pa) capped to %.1f Pa.\n', ...
            max(abs(fluid.p)), pCap);
    end
    if isfield(fluid, 'tauE') && ~isempty(fluid.tauE) && max(abs(fluid.tauE)) > tauCap
        fprintf('    [Fluid Cap Active] Peak shear stress (%.1f Pa) capped to %.1f Pa.\n', ...
            max(abs(fluid.tauE)), tauCap);
    end
    if isfield(fluid, 'tau_wall_E') && ~isempty(fluid.tau_wall_E) && max(abs(fluid.tau_wall_E)) > tauCap
        fprintf('    [Fluid Cap Active] Peak shear stress (%.1f Pa) capped to %.1f Pa.\n', ...
            max(abs(fluid.tau_wall_E)), tauCap);
    end

    % 1. Cap 1D / 2D Pressure Distributions
    if isfield(fluid, 'p') && ~isempty(fluid.p)
        fluid.p = min(max(fluid.p, -pCap), pCap);
    end
    if isfield(fluid, 'p_2d') && ~isempty(fluid.p_2d)
        fluid.p_2d = min(max(fluid.p_2d, -pCap), pCap);
    end

    % 2. Cap Wall Shear Stress Distributions
    if isfield(fluid, 'tau_wall_E') && ~isempty(fluid.tau_wall_E)
        fluid.tau_wall_E = min(max(fluid.tau_wall_E, -tauCap), tauCap);
    end
    if isfield(fluid, 'tau_wall_L') && ~isempty(fluid.tau_wall_L)
        fluid.tau_wall_L = min(max(fluid.tau_wall_L, -tauCap), tauCap);
    end

    % 3. Cap Explicit Interface Traction Vectors if constructed
    if isfield(fluid, 'TractionE_normal') && ~isempty(fluid.TractionE_normal)
        fluid.TractionE_normal = min(max(fluid.TractionE_normal, -pCap), pCap);
    end
    if isfield(fluid, 'TractionE_tangent') && ~isempty(fluid.TractionE_tangent)
        fluid.TractionE_tangent = min(max(fluid.TractionE_tangent, -tauCap), tauCap);
    end
    if isfield(fluid, 'TractionL_normal') && ~isempty(fluid.TractionL_normal)
        fluid.TractionL_normal = min(max(fluid.TractionL_normal, -pCap), pCap);
    end
    if isfield(fluid, 'TractionL_tangent') && ~isempty(fluid.TractionL_tangent)
        fluid.TractionL_tangent = min(max(fluid.TractionL_tangent, -tauCap), tauCap);
    end
end
% =========================================================================

end