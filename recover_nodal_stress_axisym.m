function stress = recover_nodal_stress_axisym(mesh, u, par)
% Recover nodal Cauchy stresses by averaging element-center stresses
% from all connected elements.

    nnode = size(mesh.nodes,1);

    sigma_rr_sum = zeros(nnode,1);
    sigma_zz_sum = zeros(nnode,1);
    sigma_rz_sum = zeros(nnode,1);
    sigma_tt_sum = zeros(nnode,1);
    vm_sum       = zeros(nnode,1);
    count        = zeros(nnode,1);

    % Use element-center stress as representative element stress
    xi = 0.0;
    eta = 0.0;
    [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        [~, dNdX, ~] = jacobian_2d(Xe, dNdxi);
        sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par);

        srr = sigma(1,1);
        stt = sigma(2,2);
        szz = sigma(3,3);
        srz = sigma(1,3);

        % axisymmetric 3D von Mises stress
        svm = sqrt(0.5*((srr-stt)^2 + (stt-szz)^2 + (szz-srr)^2 + 6*srz^2));

        for a = 1:4
            node = conn(a);
            sigma_rr_sum(node) = sigma_rr_sum(node) + srr;
            sigma_tt_sum(node) = sigma_tt_sum(node) + stt;
            sigma_zz_sum(node) = sigma_zz_sum(node) + szz;
            sigma_rz_sum(node) = sigma_rz_sum(node) + srz;
            vm_sum(node)       = vm_sum(node) + svm;
            count(node)        = count(node) + 1;
        end
    end

    stress = struct();
    stress.sigma_rr = sigma_rr_sum ./ max(count,1);
    stress.sigma_tt = sigma_tt_sum ./ max(count,1);
    stress.sigma_zz = sigma_zz_sum ./ max(count,1);
    stress.sigma_rz = sigma_rz_sum ./ max(count,1);
    stress.vonMises = vm_sum ./ max(count,1);
end

function sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par)
% Returns 3x3 Cauchy stress tensor at one quadrature point.

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    Rg = N * Rnod;
    rg = N * rnod;

    if Rg <= 0 || rg <= 0
        error('Non-positive radius encountered while post-processing stress.');
    end

    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;

    F = [drdR,   0,    drdZ;
           0,   rg/Rg, 0;
         dzdR,   0,    dzdZ];

    J = det(F);
    if J <= 0
        error('Negative or zero J encountered while post-processing stress.');
    end

    B = F * F.';
    I = eye(3);

    % Standard hyperelastic J^(-2/3) scaling (was incorrectly J^(-5/3))
    sigma = par.Ge * J^(-2/3) * ( B - (trace(B)/3)*I ) ...
          + par.Ke * (J - 1) * I;
end

