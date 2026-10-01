% CHANGES TO LOOK FOR IN THIS FILE (mainly to safeguard):
% - Lines 12-16 (Issue 1 & 3): Enforced explicit column-vectoring z(:) and interfaceNodes(:) 
%   to prevent 1xN row vector extraction from mesh.nodes and row-vector output propagation.
% - Lines 18-28 (Issue 2): Added time-step fallback guard for missing, zero, or non-finite 
%   par.dt to prevent division-by-zero (Inf/NaN) in axial wall velocity (Uw) calculations.
% - Lines 34-36 (Issue 1): Explicitly vectorized returned outputs delta and Uw to N x 1.

function [delta, Uw] = monolithic_interface_kinematics_value_only( ...
    mesh, uNew, uOld, interfaceNodes, z, par)
% MONOLITHIC_INTERFACE_KINEMATICS_VALUE_ONLY
% Value-only wrapper around monolithic_interface_kinematics.
% Maps solid interface displacements (uNew, uOld) onto fluid grid coordinates (z)
% to return gap clearances (delta) and axial wall velocities (Uw).

    if isempty(interfaceNodes) || isempty(z)
        delta = [];
        Uw = [];
        return;
    end

    % FIX 1 & 3: Explicitly sanitize input vector shapes to column vectors (N x 1)
    z_col = z(:);
    interfaceNodes_col = interfaceNodes(:);

    % FIX 2: Guard against missing, zero, or non-finite time step (dt)
    parLocal = par;
    zeroUw = false;
    if ~isfield(parLocal, 'dt') || ~isfinite(parLocal.dt) || parLocal.dt <= 0
        % If dt is invalid or zero (e.g., initial equilibrium state), set a safe dummy dt 
        % to prevent division-by-zero, then zero out Uw post-computation.
        parLocal.dt = 1.0;
        zeroUw = true;
    end

    % OLD BUGGY LINE:
    % [delta, Uw] = monolithic_interface_kinematics(mesh, uNew, uOld, interfaceNodes, z, par);

    [delta, Uw] = monolithic_interface_kinematics( ...
        mesh, uNew, uOld, interfaceNodes_col, z_col, parLocal);

    if zeroUw
        Uw = zeros(size(z_col));
    end

    % Ensure output fields maintain strict N x 1 column vector dimensions
    delta = delta(:);
    Uw = Uw(:);
end