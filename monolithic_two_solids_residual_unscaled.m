function [RE, RL, RF] = monolithic_two_solids_residual_unscaled( ...
    uE, uL, p, old, meshE, interfaceE, baseE, ... %#ok<INUSD>
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL) %#ok<INUSD>

% =========================================================================
% MANDATORY PAR / PARL SYNCHRONIZATION OVERRIDE (VERY TOP OF FILE)
% Ensures active translocation flags are never lost or overwritten
% =========================================================================
if isstruct(par) && isfield(par, 'useActiveTranslocation') && par.useActiveTranslocation
    parL.useActiveTranslocation  = true;
    if isfield(par, 'fz_active_translocation')
        parL.fz_active_translocation = par.fz_active_translocation;
    end
    if isfield(par, 'z_pore_bottom')
        parL.z_pore_bottom = par.z_pore_bottom;
    end
    if isfield(par, 'sigma_head_decay')
        parL.sigma_head_decay = par.sigma_head_decay;
    end
end
% =========================================================================

ndofE = size(meshE.nodes,1) * 2;
ndofL = size(meshL.nodes,1) * 2;
N = numel(z);

p(1) = par.pIn;
p(end) = par.pOut;

[deltaE, UwE] = monolithic_interface_kinematics_value_only( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, UwL] = monolithic_interface_kinematics_value_only( ...
    meshL, uL, old.uL, interfaceL, z, par);

gap = deltaE - deltaL;
if any(gap <= par.minGap)
    error('Gap violates minGap.');
end

[Q, ~, ~, tauL, tauE, ~, ~] = ...
    local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

[pLoadE, pLoadL] = global2d_pressure_traction_loads( ...
    z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);

% Endothelium residual. Fluid-on-endothelium traction: +p normal, -tauE tangent.
trE.normal = pLoadE;
trE.tangent = -tauE;
if use_exact_interface_in_monolithic(par)
    trE.z = z;
end

FextE = zeros(ndofE,1);
[FextE, ~] = apply_interface_traction(meshE, uE, FextE, interfaceE, trE);

FintE = assemble_finite_def_internal_force_only(meshE, uE, par) + ...
    assemble_axisym_kelvin_voigt_viscous_force_only(meshE, uE, old.uE, par);

REfull = FintE - FextE;
RE = REfull(freeE);

% Leukocyte residual. Fluid-on-leukocyte traction: -p normal, +tauL tangent.
trL.normal = -pLoadL;
trL.tangent = tauL;
if use_exact_interface_in_monolithic(par)
    trL.z = z;
end
trL = apply_leukocyte_traction_support(trL, par);

FextL = zeros(ndofL,1);
[FextL, ~] = apply_interface_traction(meshL, uL, FextL, interfaceL, trL);

FintL = assemble_finite_def_internal_force_only(meshL, uL, parL) + ...
    assemble_axisym_kelvin_voigt_viscous_force_only(meshL, uL, old.uL, parL);

% =========================================================================
% LATE-STAGE TEM ACTIVE TRANSLOCATION FORCE ASSEMBLY
% =========================================================================
% =========================================================================
% LATE-STAGE TEM ACTIVE TRANSLOCATION FORCE ASSEMBLY (CORRECTED DOF MAPPING)
% =========================================================================
f_active_global = zeros(ndofL, 1);

if isfield(parL, 'useActiveTranslocation') && parL.useActiveTranslocation
    f0      = parL.fz_active_translocation;
    z_bot   = parL.z_pore_bottom;
    sigma_z = parL.sigma_head_decay;

    % Axisymmetric FE mesh: Column 1 = r (radial), Column 2 = z (axial)
    r_nodes = meshL.nodes(:, 1);
    z_nodes = meshL.nodes(:, 2);

    % Displacement vectors: Odd DOFs = u_r, Even DOFs = u_z
    u_r = uL(1:2:end);
    u_z = uL(2:2:end);

    rL_deformed = r_nodes + u_r;
    zL_deformed = z_nodes + u_z;

    for e = 1:size(meshL.conn, 1)
        elem_nodes = meshL.conn(e, :);

        r_elem = rL_deformed(elem_nodes);
        z_elem = zL_deformed(elem_nodes);

        r_c = mean(r_elem);
        z_c = mean(z_elem);

        % Top-down translocation spatial mask:
        % Applies active propulsion force on cell body at or past z_bot (z_c <= z_bot)
        if z_c <= z_bot
            psi_z = 1.0;
        else
            psi_z = exp(-((z_c - z_bot)^2) / (2 * sigma_z^2));
        end

        area_e = polyarea(r_elem, z_elem);
        vol_e  = 2 * pi * max(r_c, 1e-9) * area_e;

        fe_z = (f0 * psi_z * vol_e) / numel(elem_nodes);

        for k = 1:numel(elem_nodes)
            node_idx = elem_nodes(k);
            % Assign propulsion force strictly to even DOFs (axial z-displacements)
            z_dof = 2 * node_idx;
            f_active_global(z_dof) = f_active_global(z_dof) + fe_z;
        end
    end
end

% =========================================================================
% UNCONSTRAINED RESIDUAL SUBTRACTION
% =========================================================================
RLfull = FintL - FextL - f_active_global;
RL     = RLfull(freeL);

% Fluid mass residual using both current interfaces.
re    = deltaE(:);
rl    = deltaL(:);
reOld = old.deltaE(:);
rlOld = old.deltaL(:);

A    = 0.5 * (re.^2    - rl.^2);
Aold = 0.5 * (reOld.^2 - rlOld.^2);

Ssrc = -par.SsrcFactor * (A - Aold) / par.dt;

dzControl = global_1d_control_lengths(z);
RF = zeros(N-2,1);
for i = 2:N-1
    RF(i-1) = (Q(i) - Q(i-1))/dzControl(i) - Ssrc(i);
end
end