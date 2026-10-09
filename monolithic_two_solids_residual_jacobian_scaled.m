function [R, J] = monolithic_two_solids_residual_jacobian_scaled( ...
    y, old, meshE, interfaceE, baseE, ...
    meshL, interfaceL, baseL, parL, ...
    z, par, freeE, fixE, valsE, freeL, fixL, valsL, ...
    JuE, JuL, Jp, solidTargetE, solidTargetL, fluidTarget)
% Semi-analytical two-solid residual/Jacobian with active translocation support.
%
% Unknown vector:
%   y = [uE_free/JuE; uL_free/JuL; p_internal/Jp]

N = numel(z);
nE = numel(freeE);
nL = numel(freeL);

try
    [uE, uL, p] = unpack_two_solid_y( ...
        y, old, freeE, fixE, valsE, freeL, fixL, valsL, JuE, JuL, Jp, par);
% Capture single output struct q from solid_geometry_quality
    qE = solid_geometry_quality(meshE, uE, 'endothelium check');
    qL = solid_geometry_quality(meshL, uL, 'leukocyte check');
    % Safe evaluation of element Jacobian minimum
    if isfield(qE, 'JEmin') && qE.JEmin <= 0.05
        error('ElementInversion:JEminTooSmall', 'Endothelium JEmin = %.3e <= 0.05', qE.JEmin);
    end
    if isfield(qL, 'JEmin') && qL.JEmin <= 0.05
        error('ElementInversion:JEminTooSmall', 'Leukocyte JEmin = %.3e <= 0.05', qL.JEmin);
    end

    % Optional: keep formal assertion wrapper if it checks additional flags (e.g. max aspect ratio)
    assert_solid_geometry_ok(qE, par);
    assert_solid_geometry_ok(qL, par);
    [RE, RL, RF, JEE, JEL, JEp, JLE, JLL, JLp, JFE, JFL, JFp] = ...
        monolithic_two_solids_residual_jacobian_unscaled( ...
        uE, uL, p, old, meshE, interfaceE, baseE, ...
        meshL, interfaceL, baseL, parL, z, par, freeE, freeL);

    R = [
        RE / solidTargetE
        RL / solidTargetL
        RF / fluidTarget
        ];

    if nargout > 1
        J = [
            (JEE*JuE)/solidTargetE, (JEL*JuL)/solidTargetE, (JEp*Jp)/solidTargetE
            (JLE*JuE)/solidTargetL, (JLL*JuL)/solidTargetL, (JLp*Jp)/solidTargetL
            (JFE*JuE)/fluidTarget,  (JFL*JuL)/fluidTarget,  (JFp*Jp)/fluidTarget
            ];
    end

    if any(~isfinite(R))
        R = 1e12 * ones(nE+nL+N-2,1);
        if nargout > 1
            J = speye(numel(R), numel(y));
        end
    end
catch ME
    if contains(ME.message, 'Gap violates minGap') || ...
            contains(ME.message, 'Negative or zero J') || ...
            contains(ME.message, 'Non-positive radius') || ...
            contains(ME.message, 'Solid geometry guard failed') || ...
            contains(ME.message, 'Element inverted')
        R = 1e12 * ones(nE+nL+N-2,1);
        if nargout > 1
            J = speye(numel(R), numel(y));
        end
    else
        rethrow(ME);
    end
end
end

function [RE, RL, RF, JEE, JEL, JEp, JLE, JLL, JLp, JFE, JFL, JFp] = ...
    monolithic_two_solids_residual_jacobian_unscaled( ...
    uE, uL, p, old, meshE, interfaceE, baseE, ... %#ok<INUSD>
    meshL, interfaceL, baseL, parL, z, par, freeE, freeL) %#ok<INUSD>

% Mandatory synchronization override
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


ndofE = size(meshE.nodes,1) * 2;
ndofL = size(meshL.nodes,1) * 2;
N = numel(z);

p(1) = par.pIn;
p(end) = par.pOut;

[deltaE, UwE, HrE, HUE] = monolithic_interface_kinematics( ...
    meshE, uE, old.uE, interfaceE, z, par);
[deltaL, UwL, HrL, HUL] = monolithic_interface_kinematics( ...
    meshL, uL, old.uL, interfaceL, z, par);

gap = deltaE - deltaL;
if any(gap <= par.minGap)
    error('Gap violates minGap.');
end

[Q, ~, ~, tauL, tauE, ~, ~] = ...
    local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

if isfield(par, 'useTwoSolidFullAnalyticalJacobian') && par.useTwoSolidFullAnalyticalJacobian
    fs = local_flux_shear_sensitivities_analytical_full(z, p, deltaL, deltaE, UwL, UwE, par);
else
    fs = local_flux_shear_sensitivities_fd_full(z, p, deltaL, deltaE, UwL, UwE, par, Q, tauL, tauE);
end

[pLoadE, pLoadL, global2DJacE, global2DJacL] = ...
    global2d_pressure_traction_loads( ...
    z, p, deltaE, deltaL, meshE, uE, meshL, uL, par);

% ------------------------------------------------------------
% Endothelium residual and same-solid tangent.
% ------------------------------------------------------------
trE.normal = pLoadE;
trE.tangent = -tauE;
if use_exact_interface_in_monolithic(par)
    trE.z = z;
end

FextE = zeros(ndofE,1);
[FextE, KextE, BnE, BtE] = ...
    apply_interface_traction_sensitivity(meshE, uE, FextE, interfaceE, trE);

[FintE, KintE] = assemble_finite_def_axisym(meshE, uE, par);
[FviscE, KviscE] = assemble_axisym_kelvin_voigt_viscous(meshE, uE, old.uE, par);
FintE = FintE + FviscE;
KintE = KintE + KviscE;

REfull = FintE - FextE;
RE = REfull(freeE);

% ------------------------------------------------------------
% Leukocyte residual and same-solid tangent.
% ------------------------------------------------------------
trL.normal = -pLoadL;
trL.tangent = tauL;
if use_exact_interface_in_monolithic(par)
    trL.z = z;
end
trL = apply_leukocyte_traction_support(trL, par);

FextL = zeros(ndofL,1);
[FextL, KextL, BnL, BtL] = ...
    apply_interface_traction_sensitivity(meshL, uL, FextL, interfaceL, trL);

[FintL, KintL] = assemble_finite_def_axisym(meshL, uL, parL);
[FviscL, KviscL] = assemble_axisym_kelvin_voigt_viscous(meshL, uL, old.uL, parL);
FintL = FintL + FviscL;
KintL = KintL + KviscL;

% =========================================================================
% LATE-STAGE TEM ACTIVE TRANSLOCATION FORCE & TANGENT ASSEMBLY
% =========================================================================
f_active_global = zeros(ndofL, 1);
K_active_global = sparse(ndofL, ndofL);

if isfield(parL, 'useActiveTranslocation') && parL.useActiveTranslocation
    f0      = parL.fz_active_translocation;
    z_bot   = parL.z_pore_bottom;
    sigma_z = parL.sigma_head_decay;
    
    r_nodes = meshL.nodes(:, 1);
    z_nodes = meshL.nodes(:, 2);
    
    % Deformed coordinates from current trial step (for element centroids)
    u_r = uL(1:2:end);
    u_z = uL(2:2:end);
    rL_deformed = r_nodes + u_r;
    zL_deformed = z_nodes + u_z;
    
    % -------------------------------------------------------------------------
    % FIX 1A-Mod: Freeze z_head using old.uL (converged displacement at t^n)
    % Prevents off-diagonal Jacobian ill-conditioning and line-search stalls
    % -------------------------------------------------------------------------
    if isstruct(old) && isfield(old, 'uL') && ~isempty(old.uL)
        u_z_old = old.uL(2:2:end);
    else
        u_z_old = u_z;
    end
    zL_deformed_old = z_nodes + u_z_old;
    if f0>=0
    z_head_frozen    = max(zL_deformed_old);
    else 
        z_head_frozen    = min(zL_deformed_old);
    end
    
    ramp_length = 0.1e-6;                             % [m] Transition length scale (100 nm)
    if f0>=0
    pore_depth  = -(z_bot - z_head_frozen);             % Positive once head penetrates past z_bot
    else
        pore_depth  = z_bot - z_head_frozen;             % Positive once head penetrates past z_bot
    end
    act_factor  = 0.5 * (1 + tanh(pore_depth / ramp_length));

    % [1008 OPT] Set 4: the tangent entries are collected as triplets and assembled with one
    % sparse() call after the loop, instead of 16 indexed updates of the sparse matrix per element.
    % Same per-element values, same order, same force-vector assembly (bit-identical, unit test).
    nTrip1008 = size(meshL.conn, 1) * size(meshL.conn, 2)^2;
    Ki1008 = zeros(nTrip1008, 1); Kj1008 = zeros(nTrip1008, 1); Kv1008 = zeros(nTrip1008, 1); p1008 = 0;

    for e = 1:size(meshL.conn, 1)
        elem_nodes = meshL.conn(e, :);
        n_elem     = numel(elem_nodes);
        
        r_elem = rL_deformed(elem_nodes);
        z_elem = zL_deformed(elem_nodes);
        
        r_c = mean(r_elem);
        z_c = mean(z_elem);
        
        % Spatial decay mask with respect to centroid z_c
        if z_c <= z_bot
            psi_z   = 1.0;
            dpsi_dz = 0.0;
        else
            psi_z   = exp(-((z_c - z_bot)^2) / (2 * sigma_z^2));
            dpsi_dz = -((z_c - z_bot) / (sigma_z^2)) * psi_z;
        end
        
        area_e = polyarea(r_elem, z_elem);
        vol_e  = 2 * pi * max(r_c, 1e-9) * area_e;
        
        % Scale active force by the frozen C1 spatial activation factor
        fe_z = (f0 * act_factor * psi_z * vol_e) / n_elem;
        
        % Active force tangent stiffness local contribution d(fe_z)/d(u_z_local)
        dfe_du_z = (f0 * act_factor * vol_e * dpsi_dz) / (n_elem^2);
        
        for k1 = 1:n_elem
            node_i = elem_nodes(k1);
            z_dof_i = 2 * node_i;
            
            % Force vector assembly
            f_active_global(z_dof_i) = f_active_global(z_dof_i) + fe_z;
            
            % Tangent stiffness assembly (sparse, local element coupling only)
            for k2 = 1:n_elem
                node_j = elem_nodes(k2);
                z_dof_j = 2 * node_j;
                p1008 = p1008 + 1;                         % [1008 OPT]
                Ki1008(p1008) = z_dof_i; Kj1008(p1008) = z_dof_j; Kv1008(p1008) = dfe_du_z;
            end
        end
    end
    K_active_global = sparse(Ki1008(1:p1008), Kj1008(1:p1008), Kv1008(1:p1008), ndofL, ndofL);   % [1008 OPT]

end

        
% Subtract active force from leukocyte residual
RLfull = FintL - FextL - f_active_global;
RL     = RLfull(freeL);

% ------------------------------------------------------------
% Fluid residual.
% ------------------------------------------------------------
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

Dq = sparse(N-2, N-1);
for i = 2:N-1
    row = i-1;
    Dq(row,i)   =  1/dzControl(i);
    Dq(row,i-1) = -1/dzControl(i);
end

KFp_all = Dq * fs.dQdp;
KFre = Dq * fs.dQdre;
KFrl = Dq * fs.dQdrl;
KFUwE = Dq * fs.dQdUwE;
KFUwL = Dq * fs.dQdUwL;

for i = 2:N-1
    row = i-1;
    KFre(row,i) = KFre(row,i) + par.SsrcFactor * re(i) / par.dt;
    KFrl(row,i) = KFrl(row,i) - par.SsrcFactor * rl(i) / par.dt;
end

JFEfull = KFre*HrE + KFUwE*HUE;
JFLfull = KFrl*HrL + KFUwL*HUL;
JFE = JFEfull(:,freeE);
JFL = JFLfull(:,freeL);
JFp = KFp_all(:,2:end-1);

% ------------------------------------------------------------
% Solid Jacobian blocks.
% ------------------------------------------------------------
dtauE_duE = fs.dtauEdre*HrE + fs.dtauEdUwE*HUE;
dtauE_duL = fs.dtauEdrl*HrL + fs.dtauEdUwL*HUL;
dtauL_duE = fs.dtauLdre*HrE + fs.dtauLdUwE*HUE;
dtauL_duL = fs.dtauLdrl*HrL + fs.dtauLdUwL*HUL;

JEEfull = KintE - KextE + BtE*dtauE_duE;
JELfull =              + BtE*dtauE_duL;
JEpfull = -BnE*global2DJacE(:,2:end-1) + BtE*fs.dtauEdp(:,2:end-1);

% Leukocyte tangent: R_L = F_int - F_ext - f_active -> J_LL = K_int - K_ext - d(f_active)/du_L
JLEfull =              - BtL*dtauL_duE;
JLLfull = KintL - KextL - BtL*dtauL_duL - K_active_global;
JLpfull =  BnL*global2DJacL(:,2:end-1) - BtL*fs.dtauLdp(:,2:end-1);

JEE = JEEfull(freeE, freeE);
JEL = JELfull(freeE, freeL);
JEp = JEpfull(freeE, :);

JLE = JLEfull(freeL, freeE);
JLL = JLLfull(freeL, freeL);
JLp = JLpfull(freeL, :);
end