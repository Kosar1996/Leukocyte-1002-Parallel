function stress = recover_nodal_stress_axisym_viscoelastic(mesh, u, uOld, dt, par)
% RECOVER_NODAL_STRESS_AXISYM_VISCOELASTIC
% Recovers smoothed nodal viscoelastic stresses across axisymmetric elements.
%
% Modified to use 4-point Gauss quadrature, 2x2 Gauss-to-node extrapolation,
% and von Mises invariant calculation while maintaining the external call to
% objective_kelvin_voigt_piola.m for the viscous stress engine.

nnode = size(mesh.nodes,1);

sigma_rr_sum = zeros(nnode,1);
sigma_zz_sum = zeros(nnode,1);
sigma_rz_sum = zeros(nnode,1);
sigma_tt_sum = zeros(nnode,1);
% [NEW]: Added nodal summation array for von Mises stress invariant
vm_sum       = zeros(nnode,1);
count        = zeros(nnode,1);

% [OLD]: xi = 0.0; eta = 0.0;
% [OLD]: [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);

% [NEW]: Local node locations for 2x2 Q4 element extrapolation
xi_nodes  = [-1,  1,  1, -1];
eta_nodes = [-1, -1,  1,  1];

% [NEW]: Standard 2x2 Gauss quadrature points
g_val = 1.0 / sqrt(3);
gps = [-g_val, -g_val;
        g_val, -g_val;
        g_val,  g_val;
       -g_val,  g_val];

useVisc = isfield(par,'etaE') && isfinite(par.etaE) && par.etaE > 0;

parVisc = par;
parVisc.dt = dt;

for e = 1:mesh.nelem
    conn = mesh.conn(e,:);
    Xe   = mesh.nodes(conn,:);
    dofs = reshape([2*conn-1; 2*conn], [], 1);
    ue    = u(dofs);
    ueOld = uOld(dofs);

    % [NEW]: Struct to hold stresses evaluated across the 4 Gauss integration points
    gp_stress = struct();
    gp_stress.srr = zeros(4,1);
    gp_stress.stt = zeros(4,1);
    gp_stress.szz = zeros(4,1);
    gp_stress.srz = zeros(4,1);
    gp_stress.svm = zeros(4,1);

    % [OLD]: [~, dNdX, ~] = jacobian_2d(Xe, dNdxi);
    % [OLD]: F = local_F(Xe, ue, N, dNdX);
    % [OLD]: J = det(F);
    % [OLD]: if J <= 0, error('Negative or zero J encountered while post-processing stress.'); end
    % [OLD]: B = F * F.';
    % [OLD]: I3 = eye(3);
    % [OLD]: sigmaElastic = par.Ge * J^(-2/3) * ( B - (trace(B)/3)*I3 ) + par.Ke * (J - 1) * I3;
    % [OLD]: if useVisc
    % [OLD]:     Fold = local_F(Xe, ueOld, N, dNdX);
    % [OLD]:     [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
    % [OLD]:     sigma = sigmaElastic + kvdata.Tvisc;
    % [OLD]: else
    % [OLD]:     sigma = sigmaElastic;
    % [OLD]: end
    % [OLD]: srr = sigma(1,1); stt = sigma(2,2); szz = sigma(3,3); srz = sigma(1,3);

    % [NEW]: 1. Full 4-Point Gauss Quadrature Evaluation
    for g = 1:4
        xi  = gps(g,1);
        eta = gps(g,2);
        
        [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
        [~, dNdX, ~]  = jacobian_2d(Xe, dNdxi);

        F = local_F(Xe, ue, N, dNdX);
        J = det(F);
        if J <= 0
            error('Negative or zero J encountered while post-processing stress.');
        end
        
        B = F * F.';
        I3 = eye(3);
        % Standard hyperelastic J^(-2/3) scaling
        sigmaElastic = par.Ge * J^(-2/3) * ( B - (trace(B)/3)*I3 ) + par.Ke * (J - 1) * I3;

        % External call to objective_kelvin_voigt_piola engine maintained
        if useVisc
            Fold = local_F(Xe, ueOld, N, dNdX);
            [~, kvdata] = objective_kelvin_voigt_piola(F, Fold, parVisc);
            sigma = sigmaElastic + kvdata.Tvisc;
        else
            sigma = sigmaElastic;
        end

        srr = sigma(1,1);
        stt = sigma(2,2);
        szz = sigma(3,3);
        srz = sigma(1,3);

        % [NEW]: 2. von Mises Stress Invariant Calculation
        svm = sqrt(0.5 * ((srr - stt)^2 + (stt - szz)^2 + (szz - srr)^2 + 6 * srz^2));

        gp_stress.srr(g) = srr;
        gp_stress.stt(g) = stt;
        gp_stress.szz(g) = szz;
        gp_stress.srz(g) = srz;
        gp_stress.svm(g) = svm;
    end

    % [OLD]: Direct single-point element averaging to nodes:
    % [OLD]: for a = 1:4
    % [OLD]:     node = conn(a);
    % [OLD]:     sigma_rr_sum(node) = sigma_rr_sum(node) + srr;
    % [OLD]:     sigma_tt_sum(node) = sigma_tt_sum(node) + stt;
    % [OLD]:     sigma_zz_sum(node) = sigma_zz_sum(node) + szz;
    % [OLD]:     sigma_rz_sum(node) = sigma_rz_sum(node) + srz;
    % [OLD]:     count(node)        = count(node) + 1;
    % [OLD]: end

    % [NEW]: 3. Extrapolation from 2x2 Gauss Points to Element Corner Nodes
    a_sqrt3 = sqrt(3);
    for a_idx = 1:4
        xi_n = xi_nodes(a_idx);
        eta_n = eta_nodes(a_idx);

        % Extrapolation shape functions evaluated at nodal positions
        E = 0.25 * [ (1 - a_sqrt3*xi_n)*(1 - a_sqrt3*eta_n), ...
                     (1 + a_sqrt3*xi_n)*(1 - a_sqrt3*eta_n), ...
                     (1 + a_sqrt3*xi_n)*(1 + a_sqrt3*eta_n), ...
                     (1 - a_sqrt3*xi_n)*(1 + a_sqrt3*eta_n) ];

        node = conn(a_idx);
        sigma_rr_sum(node) = sigma_rr_sum(node) + E * gp_stress.srr;
        sigma_tt_sum(node) = sigma_tt_sum(node) + E * gp_stress.stt;
        sigma_zz_sum(node) = sigma_zz_sum(node) + E * gp_stress.szz;
        sigma_rz_sum(node) = sigma_rz_sum(node) + E * gp_stress.srz;
        vm_sum(node)       = vm_sum(node)       + E * gp_stress.svm;
        count(node)        = count(node) + 1;
    end
end

% stress = struct();
% stress.sigma_rr = sigma_rr_sum ./ max(count,1);
% stress.sigma_tt = sigma_tt_sum ./ max(count,1);
% stress.sigma_zz = sigma_zz_sum ./ max(count,1);
% stress.sigma_rz = sigma_rz_sum ./ max(count,1);
% % [NEW]: Added von Mises field to output struct
% stress.vonMises = vm_sum       ./ max(count,1);

% ========================= NEW CODE ==============================
s_rr = sigma_rr_sum ./ max(count,1);
s_tt = sigma_tt_sum ./ max(count,1);
s_zz = sigma_zz_sum ./ max(count,1);
s_rz = sigma_rz_sum ./ max(count,1);
s_vm = vm_sum       ./ max(count,1);

% Enforce Axisymmetric Boundary Regularity exactly at r = 0
axis_nodes = find(mesh.nodes(:,1) < 1e-10);
if ~isempty(axis_nodes)
    % 1. Zero shear stress at symmetry axis
    s_rz(axis_nodes) = 0.0;
    
    % 2. Enforce equal radial and hoop stress (sigma_rr = sigma_tt)
    s_avg = 0.5 * (s_rr(axis_nodes) + s_tt(axis_nodes));
    s_rr(axis_nodes) = s_avg;
    s_tt(axis_nodes) = s_avg;
    
    % 3. Re-evaluate von Mises on axis using regularized stresses
    s_vm(axis_nodes) = sqrt(0.5 * ((s_rr(axis_nodes) - s_tt(axis_nodes)).^2 + ...
                                  (s_tt(axis_nodes) - s_zz(axis_nodes)).^2 + ...
                                  (s_zz(axis_nodes) - s_rr(axis_nodes)).^2 + ...
                                  6 * s_rz(axis_nodes).^2));
end

stress = struct();
stress.sigma_rr = s_rr;
stress.sigma_tt = s_tt;
stress.sigma_zz = s_zz;
stress.sigma_rz = s_rz;
stress.vonMises  = s_vm;
% =================================================================

end

function F = local_F(Xe, ue, N, dNdX)
Rnod = Xe(:,1);
Znod = Xe(:,2);
rnod = Rnod + ue(1:2:end);
znod = Znod + ue(2:2:end);

Rg = N * Rnod;
rg = N * rnod;
% if Rg <= 0 || rg <= 0
%     error('Non-positive radius encountered while post-processing stress.');
% end
% 
% drdR = dNdX(:,1).' * rnod;
% drdZ = dNdX(:,2).' * rnod;
% dzdR = dNdX(:,1).' * znod;
% dzdZ = dNdX(:,2).' * znod;
% 
% F = [drdR,   0,    drdZ;
%        0,   rg/Rg, 0;
%      dzdR,   0,    dzdZ];

% ========================= NEW CODE ==============================
% Centerline Regularization & L'Hopital Safeguards
epsR = 1e-14;
Rg_eff = max(Rg, epsR);
rg_eff = max(rg, epsR);

drdR = dNdX(:,1).' * rnod;
drdZ = dNdX(:,2).' * rnod;
dzdR = dNdX(:,1).' * znod;
dzdZ = dNdX(:,2).' * znod;

% L'Hopital Limit for Hoop Stretch on axis
if Rg < 1e-10
    F22 = drdR;
else
    F22 = rg_eff / Rg_eff;
end

F = [drdR,   0,    drdZ;
       0,   F22,   0;
     dzdR,   0,    dzdZ];
% =================================================================
end