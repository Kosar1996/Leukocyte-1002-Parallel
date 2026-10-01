clc;clear;close all;
% =========================
par = struct();

par.RLin     = 0e-6;       % leukocyte inner radius [m]
par.RLout    = 4e-6;       % leukocyte outer radius [m]

par.Rc      = 2e-6;

par.zMinFluid = 0e-6;         % coupled fluid lower axial bound [m]
par.zMaxFluid = 4e-6;         % coupled fluid upper axial bound [m]
par.LzFluid   = par.zMaxFluid - par.zMinFluid;

par.zMin    = -2e-6;          % leukocyte lower axial bound [m] (FIXED: typo -2-6)
par.zMax    = 6e-6;           % leukocyte upper axial bound [m]
par.Lz      = par.zMax - par.zMin;

par.NzFluid = 61;             % coupled fluid grid points
par.NrL     = 30;             % radial nodes, leukocyte
par.NzSolid = 121;

par.EL  = 200;       % leukocyte Young's modulus [Pa]
par.nuL = 0.46;      % leukocyte Poisson ratio

par.Gl = par.EL/(2*(1+par.nuL));
par.Kl = par.EL/(3*(1-2*par.nuL));
par.useViscoelasticLeukocyte = true;

par.pIn     = 0;            % inlet pressure [Pa]
par.pOut    = 0;            % outlet pressure [Pa]

par.supportL = 'roller';    % axis, clamped, or roller

par.newtonMaxItSolid = 100;
par.lineSearchMax    = 50;
par.newtonTolSolid   = 1e-8;
par.nLoadStepsSolid  = 20;

% Building the finite-element meshes for the solid domain
meshL = build_rounded_leukocyte_mesh(par);
interfaceL = find_interface_nodes(meshL, 'outer');  % facing fluid gap
baseL      = find_interface_nodes(meshL, 'inner');  % central axis/anchored side

zSolid = linspace(par.zMin, par.zMax, par.NzSolid).';
zFluid = linspace(par.zMinFluid, par.zMaxFluid, par.NzFluid).';

par.P0l = 600;                 % Pa
tau0 = 0;                      % Pa

zLoad = zSolid;   % use the leukocyte axial grid for the prestress load
zLoadMin = 0e-6;
zLoadMax = 4e-6;

%% Step-like pressure profile
edgeWidth = 0.5e-6;
leftEdge  = 0.5 * (1 + tanh((zLoad - zLoadMin) / edgeWidth));
rightEdge = 0.5 * (1 - tanh((zLoad - zLoadMax) / edgeWidth));
pressureProfile = leftEdge .* rightEdge;

%% Cos-like pressure profile
zMid = 0.5*(zLoadMin + zLoadMax);
halfWidth = 0.5*(zLoadMax - zLoadMin);
xi = (zLoad - zMid)/halfWidth;
pressureProfile = 0.5*(1 + cos(pi*xi));
pressureProfile(abs(xi) > 1) = 0;

trL.normal  = -par.P0l * pressureProfile;
trL.tangent = tau0 * ones(size(zLoad));   % uniform shear

% =========================

state = initial_state(zSolid, par, meshL);
uLnew = solve_finite_def_leukocyte(meshL, state.uE, trL, interfaceL, baseL, par);

%Postprocess interface shape
deltal = extract_interface_radius(meshL, uLnew, interfaceL, zLoad);
uzL    = extract_interface_axial_displacement(meshL, uLnew, interfaceL, zLoad);

figure;
hold on;
plot_deformed_axisym_mesh(meshL, uLnew, 1);   % scale = 1 means exact deformed position
xlabel('r [\mum]');
ylabel('z [\mum]');
axis equal;
set(gca,'FontSize',14);
box on;
xlim([0 15]);

% Plot radial position of outer surface
figure;
plot(zLoad*1e6, deltal*1e6, 'LineWidth', 2);
xlabel('z [\mum]');
ylabel('deformed outer radius [\mum]');
title('Leukocyte outer surface');

surfL = compute_surface_traction(meshL, uLnew, interfaceL, 'outer', par);

% Recover nodal stresses using 4-point bilinear extrapolation
stressL = recover_nodal_stress_axisym(meshL, uLnew, par);

% Plot stress contours
plot_nodal_stress_contour(meshL, uLnew, stressL.sigma_rr, ...
    '\sigma_{rr} [Pa]', 'Leukocyte nodal stress contour: \sigma_{rr}');

plot_nodal_stress_contour(meshL, uLnew, stressL.sigma_zz, ...
    '\sigma_{zz} [Pa]', 'Leukocyte nodal stress contour: \sigma_{zz}');

plot_nodal_stress_contour(meshL, uLnew, stressL.sigma_rz, ...
    '\sigma_{rz} [Pa]', 'Leukocyte nodal stress contour: \sigma_{rz}');

plot_nodal_stress_contour(meshL, uLnew, stressL.vonMises, ...
    '\sigma_{vM} [Pa]', 'Leukocyte nodal stress contour: von Mises');

plot_nodal_stress_contour(meshL, uLnew, stressL.sigma_tt, ...
    '\sigma_{\theta\theta} [Pa]', ...
    'Leukocyte nodal stress contour: \sigma_{\theta\theta}');

uL_pre = uLnew;

deltaL_pre = extract_interface_radius(meshL, uL_pre, interfaceL, zFluid);
z = zFluid;

save('prestress2.mat', 'par', 'uL_pre', 'deltaL_pre', 'z', 'meshL', 'interfaceL', 'baseL', 'zSolid', 'zFluid');


%% ===================== HELPER & SOLVER FUNCTIONS =====================

function state = initial_state(z, par, meshE)
    state = struct();
    state.UwL = zeros(size(z));
    state.UwE = zeros(size(z));

    state.deltaL = rounded_outer_profile(z, par.RLout, par.zMin, par.Lz, par.Rc);
    state.deltaE = zeros(size(z));

    state.p = linspace(par.pIn, par.pOut, numel(z)).';

    state.uE = zeros(size(meshE.nodes,1)*2,1);
end


function ids = find_interface_nodes(mesh, whichSide)
    zvals = unique(mesh.nodes(:,2));
    nz = numel(zvals);
    nnode = size(mesh.nodes,1);
    nr = nnode / nz;

    if abs(nr - round(nr)) > 1e-12
        error('Cannot infer structured mesh dimensions.');
    end

    nr = round(nr);

    switch lower(whichSide)
        case 'inner'
            i = 1;
        case 'outer'
            i = nr;
        otherwise
            error('unknown side');
    end

    j = (1:nz).';
    ids = sub2ind([nz,nr], j, i*ones(nz,1));
end


function uzInt = extract_interface_axial_displacement(mesh, u, interfaceNodes, zq)
    zn = mesh.nodes(interfaceNodes,2);
    uz = u(2*interfaceNodes);
    [zs, idx] = sort(zn);
    uzs = uz(idx);
    uzInt = interp1(zs, uzs, zq, 'linear', 'extrap');
end


function [fixDofs, fixVals] = solid_support_conditions(baseNodes, supportType)
    switch lower(supportType)
        case 'axis'
            fixDofs = sort(2*baseNodes(:)-1);
            n0 = baseNodes(round(end/2));
            fixDofs = sort([fixDofs; 2*n0]);
            fixVals = zeros(size(fixDofs));

        case 'clamped'
            fixDofs = sort([2*baseNodes(:)-1; 2*baseNodes(:)]);
            fixVals = zeros(size(fixDofs));

        case 'roller'
            fixDofs = sort(2*baseNodes(:)-1);
            n0 = baseNodes(round(end/2));
            fixDofs = sort([fixDofs; 2*n0]);
            fixVals = zeros(size(fixDofs));

        otherwise
            error('Unknown support type: %s', supportType);
    end
end


function delta = extract_interface_radius(mesh, u, interfaceNodes, zq)
    zn = mesh.nodes(interfaceNodes,2);
    rn = mesh.nodes(interfaceNodes,1);
    un = u(2*interfaceNodes-1);

    [zs, idx] = sort(zn);
    rs = rn(idx) + un(idx);

    %fix added the change in zs to prevent aphysical overlap of delta
    %values
    uz = u(2*interfaceNodes);
    zs=zs+uz(idx);

    delta = interp1(zs, rs, zq, 'linear', 'extrap');
end


function zq = traction_to_z(traction, zNodes)
    N = numel(traction.normal);
    z0 = linspace(min(zNodes), max(zNodes), N).';
    zq.normal  = interp1(z0, traction.normal,  zNodes, 'linear', 'extrap');
    zq.tangent = interp1(z0, traction.tangent, zNodes, 'linear', 'extrap');
end


function [N, dNdxi, w] = q4_shape(xi, eta, gw)
    N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dNdxi = 0.25*[ -(1-eta), -(1-xi);
                    +(1-eta), -(1+xi);
                    +(1+eta), +(1+xi);
                    -(1+eta), +(1-xi) ];
    w = gw;
end


function [Jm, dNdx, detJ] = jacobian_2d(Xe, dNdxi)
    Jm = Xe.' * dNdxi;
    detJ = det(Jm);
    if detJ <= 0
        error('Non-positive element Jacobian.');
    end
    dNdx = dNdxi / Jm;
end


function uNew = solve_finite_def_leukocyte(mesh, uOld, traction, interfaceNodes, baseNodes, par)
    ndof = size(mesh.nodes,1)*2;
    u = uOld;

    [fixDofs, fixVals] = solid_support_conditions(baseNodes, par.supportL);
    free = setdiff((1:ndof).', unique(fixDofs(:)));
    u(fixDofs) = fixVals;

    nLoadSteps = par.nLoadStepsSolid;
    tractionFull = traction;

    for loadStep = 1:nLoadSteps
        scale = loadStep/nLoadSteps;
        traction.normal = scale * tractionFull.normal;
        traction.tangent = scale * tractionFull.tangent;
        fprintf('Load step %d/%d\n', loadStep, nLoadSteps);

        for it = 1:par.newtonMaxItSolid
            Fext = zeros(ndof,1);
            [Fext, Kext] = apply_interface_traction(mesh, u, Fext, interfaceNodes, traction);
        
            [Fint, Ktan] = assemble_finite_def_axisym(mesh, u, par);
        
            R = Fint - Fext;
            Ktot = Ktan - Kext;
        
            Rf  = R(free);
            Kff = Ktot(free, free);
            fprintf('   it=%02d rcond(Kff)=%.3e\n', it, rcond(full(Kff)));

            resNorm = norm(Rf, inf);
            refNorm = max(norm(Fext(free), inf), 1e-14);
            relRes = resNorm/refNorm;
            fprintf('   relRes=%.3e\n', relRes);
            if relRes < par.newtonTolSolid
                break;
            end

            du_free = -Kff \ Rf;

            fprintf('   max|du_free|=%.3e, min|du_free|=%.3e\n', max(abs(du_free)), min(abs(du_free)));

            alpha = 1.0;
            accepted = false;

            for ls = 1:par.lineSearchMax
                uTrial = u;
                uTrial(free) = uTrial(free) + alpha*du_free;
                uTrial(fixDofs) = fixVals;
        
                try
                    FextTrial = zeros(ndof,1);
                    [FextTrial, ~] = apply_interface_traction(mesh, uTrial, FextTrial, interfaceNodes, traction);
        
                    FintTrial = assemble_finite_def_internal_force_only(mesh, uTrial, par);
                    Rtrial = FintTrial - FextTrial;
        
                    if norm(Rtrial(free), inf) < resNorm
                        u = uTrial;
                        accepted = true;
                        fprintf('   line search it=%02d\n', ls);
                        break;
                    end
        
                catch ME
                    if contains(ME.message, 'Negative or zero J') || ...
                       contains(ME.message, 'Non-positive radius') || ...
                       contains(ME.message, 'Element inverted')
                    else
                        rethrow(ME);
                    end
                end
        
                alpha = 0.5 * alpha;
            end
        
            if ~accepted
                warning('Line search failed: no non-inverted residual-reducing step found.');
                uNew = u;
                return;
            end
        end

        if relRes >= par.newtonTolSolid
            warning('Finite-deformation solve hit max Newton iterations at load step %d.', loadStep);
            uNew = u;
            return;
        end
    end

    uNew = u;
end


function [Fint, K] = assemble_finite_def_axisym(mesh, u, par)
    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);
    K    = sparse(ndof, ndof);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        [fe, Ke] = finite_def_element_residual_tangent(Xe, ue, mesh, par);

        Fint(dofs) = Fint(dofs) + fe;
        K(dofs,dofs) = K(dofs,dofs) + Ke;
    end
end


function Fint = assemble_finite_def_internal_force_only(mesh, u, par)
    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        fe = finite_def_element_residual_only(Xe, ue, mesh, par);
        Fint(dofs) = Fint(dofs) + fe;
    end
end


function [fe, Ke] = finite_def_element_residual_tangent(Xe, ue, mesh, par)
    fe = zeros(8,1);
    Ke = zeros(8,8);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

    I3 = eye(3);

    for g = 1:mesh.ngp
        xi  = mesh.gp(g,1);
        eta = mesh.gp(g,2);
        w   = mesh.gw(g);

        [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
        [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

        Rg = N * Rnod;
        rg = N * rnod;

        if Rg <= 0 || rg <= 0
            error('Non-positive radius encountered in finite-deformation element.');
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
            error('Negative or zero J encountered. Element inverted.');
        end

        Finv  = inv(F);
        FinvT = Finv.';
        
        % Standard hyperelastic J^(-2/3) scaling
        J23  = J^(-2/3);
        B    = F * F.';
        trB  = trace(B);
        devB = B - (trB/3)*I3;

        % Cauchy Stress
        T = par.Gl * J23 * devB + par.Kl * (J - 1) * I3;

        % First Piola-Kirchhoff Stress
        P = J * T * FinvT;

        % Reference domain weight
        Wgp = (2*pi*Rg) * detJ0 * w;

        % ---- residual contribution ----
        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) * Wgp;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
        end

        % ---- consistent analytical tangent ----
        for alpha = 1:8
            dF = local_dF_from_dof(alpha, N, dNdX, Rg);

            dJ = J * trace(Finv * dF);

            dB = dF * F.' + F * dF.';
            trdB = trace(dB);
            dDevB = dB - (trdB/3)*I3;

            dJ23 = -(2/3) * J23 * trace(Finv * dF);

            dT = par.Gl * ( dJ23 * devB + J23 * dDevB ) ...
               + par.Kl * dJ * I3;

            dFinvT = -FinvT * dF.' * FinvT;

            dP = dJ * T * FinvT + J * dT * FinvT + J * T * dFinvT;

            for a = 1:4
                dNa_dR = dNdX(a,1);
                dNa_dZ = dNdX(a,2);
                Na     = N(a);

                Ke(2*a-1, alpha) = Ke(2*a-1, alpha) + ...
                    ( dP(1,1)*dNa_dR + dP(1,3)*dNa_dZ + dP(2,2)*(Na/Rg) ) * Wgp;

                Ke(2*a, alpha) = Ke(2*a, alpha) + ...
                    ( dP(3,1)*dNa_dR + dP(3,3)*dNa_dZ ) * Wgp;
            end
        end
    end
end


function dF = local_dF_from_dof(alpha, N, dNdX, Rg)
    aNode = ceil(alpha/2);
    isRadial = mod(alpha,2)==1;

    dF = zeros(3,3);

    if isRadial
        d_r_g   = N(aNode);
        d_drdR  = dNdX(aNode,1);
        d_drdZ  = dNdX(aNode,2);

        dF(1,1) = d_drdR;
        dF(1,3) = d_drdZ;
        dF(2,2) = d_r_g / Rg;
    else
        d_dzdR = dNdX(aNode,1);
        d_dzdZ = dNdX(aNode,2);

        dF(3,1) = d_dzdR;
        dF(3,3) = d_dzdZ;
    end
end


function fe = finite_def_element_residual_only(Xe, ue, mesh, par)
    fe = zeros(8,1);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);
    
    I3 = eye(3);

    for g = 1:mesh.ngp
        xi  = mesh.gp(g,1);
        eta = mesh.gp(g,2);
        w   = mesh.gw(g);

        [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
        [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

        Rg = N * Rnod;
        rg = N * rnod;

        if Rg <= 0 || rg <= 0
            error('Non-positive radius encountered in finite-deformation element.');
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
            error('Negative or zero J encountered. Element inverted.');
        end

        B = F * F.';

        % Standard hyperelastic Cauchy stress matching tangent implementation
        T = par.Gl * J^(-2/3) * ( B - (trace(B)/3)*I3 ) ...
          + par.Kl * (J - 1) * I3;

        P = J * T / (F.');
        Wgp = (2*pi*Rg) * detJ0 * w;

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) * Wgp;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) * Wgp;
        end
    end
end


function [F, Kext] = apply_interface_traction(mesh, u, F, interfaceNodes, traction)
    ndof = size(mesh.nodes,1) * 2;
    Kext = sparse(ndof, ndof);

    zn = mesh.nodes(interfaceNodes,2);
    [zs, idx] = sort(zn);
    interface = interfaceNodes(idx);

    zq = traction_to_z(traction, zs);
    tr_n = zq.normal;
    tr_t = zq.tangent;

    xi_gp = [-1, 1] / sqrt(3);
    w_gp  = [1, 1];

    for k = 1:numel(interface)-1
        n1 = interface(k);
        n2 = interface(k+1);

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

        r1 = mesh.nodes(n1,1);  z1 = mesh.nodes(n1,2);
        r2 = mesh.nodes(n2,1);  z2 = mesh.nodes(n2,2);

        ur1 = u(2*n1-1); uz1 = u(2*n1);
        ur2 = u(2*n2-1); uz2 = u(2*n2);

        x1 = [r1 + ur1; z1 + uz1];
        x2 = [r2 + ur2; z2 + uz2];

        dx = x2 - x1;
        L  = norm(dx);

        if L <= 1e-12
            continue;
        end

        t_hat = dx / L;
        n_hat = [ t_hat(2); -t_hat(1) ];

        tn_nodes = [tr_n(k);   tr_n(k+1)];
        tt_nodes = [tr_t(k);   tr_t(k+1)];

        fe = zeros(4,1);
        ke = zeros(4,4);

        for g = 1:2
            xi = xi_gp(g);
            wg = w_gp(g);

            N1 = 0.5 * (1 - xi);
            N2 = 0.5 * (1 + xi);

            Nline = [N1, N2];

            r_gp = N1 * x1(1) + N2 * x2(1);

            tn_gp = Nline * tn_nodes;
            tt_gp = Nline * tt_nodes;

            tvec = tn_gp * n_hat + tt_gp * t_hat;

            Nmat = [N1 0  N2 0;
                    0  N1 0  N2];

            Jline = L / 2;
            fac   = (2*pi*r_gp) * Jline * wg;

            fe = fe + (Nmat.' * tvec) * fac;

            Bdx = [-1  0  1  0;
                    0 -1  0  1];

            Br  = [N1 0 N2 0];

            I2 = eye(2);
            Ptan = I2 - (t_hat * t_hat.');
            R90 = [0 1; -1 0];

            At = (Ptan / L) * Bdx;
            An = R90 * At;

            AJ = 0.5 * (t_hat.' * Bdx);
            Afac = 2*pi * wg * ( Jline * Br + r_gp * AJ );
            Avec = tn_gp * An + tt_gp * At;

            ke = ke + (Nmat.' * Avec) * fac + (Nmat.' * tvec) * Afac;
        end

        F(dofs) = F(dofs) + fe;
        Kext(dofs,dofs) = Kext(dofs,dofs) + ke;
    end
end


function surf = compute_surface_traction(mesh, u, interfaceNodes, whichSide, par)
    zn = mesh.nodes(interfaceNodes,2);
    [~, idx] = sort(zn);
    interface = interfaceNodes(idx);

    xi_gp = [-1, 1] / sqrt(3);
    w_gp  = [1, 1];

    z_all  = [];
    r_all  = [];
    tn_all = [];
    tt_all = [];
    tr_all = [];
    tz_all = [];
    srr_all = [];
    srz_all = [];
    szz_all = [];

    Fr_total = 0;
    Fz_total = 0;
    Fn_total = 0;
    Ft_total = 0;

    for k = 1:numel(interface)-1
        n1 = interface(k);
        n2 = interface(k+1);

        e = find_boundary_element(mesh, n1, n2, whichSide);
        if isempty(e)
            continue;
        end

        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        xcurr = Xe;
        xcurr(:,1) = xcurr(:,1) + ue(1:2:end);
        xcurr(:,2) = xcurr(:,2) + ue(2:2:end);

        switch lower(whichSide)
            case 'inner'
                xi_fixed = -1;
            case 'outer'
                xi_fixed = +1;
            otherwise
                error('whichSide must be ''inner'' or ''outer''.');
        end

        for g = 1:2
            eta = xi_gp(g);
            wg  = w_gp(g);

            [N, dNdxi, ~] = q4_shape(xi_fixed, eta, 1.0);
            [~, dNdX, ~]  = jacobian_2d(Xe, dNdxi);

            rg = N * xcurr(:,1);
            zg = N * xcurr(:,2);

            dN_deta = 0.25 * [ -(1-xi_fixed);
                               -(1+xi_fixed);
                               +(1+xi_fixed);
                               +(1-xi_fixed) ];
            dx_deta = [dN_deta.' * xcurr(:,1);
                       dN_deta.' * xcurr(:,2)];

            ds_deta = norm(dx_deta);
            if ds_deta <= 0
                continue;
            end

            s_hat = dx_deta / ds_deta;
            n_hat_base = [ s_hat(2); -s_hat(1) ];

            if strcmpi(whichSide, 'inner')
                n_hat = -n_hat_base;
            else
                n_hat = n_hat_base;
            end

            sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par);

            srr = sigma(1,1);
            srz = sigma(1,3);
            szz = sigma(3,3);

            tvec = [srr, srz;
                    srz, szz] * n_hat;

            t_hat = s_hat;

            tn = dot(tvec, n_hat);
            tt = dot(tvec, t_hat);

            dA = 2*pi*rg*ds_deta*wg;

            Fr_total = Fr_total + tvec(1) * dA;
            Fz_total = Fz_total + tvec(2) * dA;
            Fn_total = Fn_total + tn      * dA;
            Ft_total = Ft_total + tt      * dA;

            z_all   = [z_all; zg];
            r_all   = [r_all; rg];
            tn_all  = [tn_all; tn];
            tt_all  = [tt_all; tt];
            tr_all  = [tr_all; tvec(1)];
            tz_all  = [tz_all; tvec(2)];
            srr_all = [srr_all; srr];
            srz_all = [srz_all; srz];
            szz_all = [szz_all; szz];
        end
    end

    [z_all, ord] = sort(z_all);

    surf = struct();
    surf.z = z_all;
    surf.r = r_all(ord);

    surf.tn = tn_all(ord);
    surf.tt = tt_all(ord);
    surf.tr = tr_all(ord);
    surf.tz = tz_all(ord);

    surf.sigma_rr = srr_all(ord);
    surf.sigma_rz = srz_all(ord);
    surf.sigma_zz = szz_all(ord);

    surf.Fr_total = Fr_total;
    surf.Fz_total = Fz_total;
    surf.Fn_total = Fn_total;
    surf.Ft_total = Ft_total;
end


function e = find_boundary_element(mesh, n1, n2, whichSide)
    e = [];

    for ii = 1:size(mesh.conn,1)
        conn = mesh.conn(ii,:);

        switch lower(whichSide)
            case 'inner'
                edge = conn([1 4]);
            case 'outer'
                edge = conn([2 3]);
            otherwise
                error('whichSide must be ''inner'' or ''outer''.');
        end

        if isequal(edge(:), [n1; n2]) || isequal(edge(:), [n2; n1])
            e = ii;
            return;
        end
    end
end


function sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par)
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
    I3 = eye(3);

    sigma = par.Gl * J^(-2/3) * ( B - (trace(B)/3)*I3 ) ...
          + par.Kl * (J - 1) * I3;
end


function stress = recover_nodal_stress_axisym(mesh, u, par)
    nnode = size(mesh.nodes,1);

    sigma_rr_sum = zeros(nnode,1);
    sigma_zz_sum = zeros(nnode,1);
    sigma_rz_sum = zeros(nnode,1);
    sigma_tt_sum = zeros(nnode,1);
    vm_sum       = zeros(nnode,1);
    count        = zeros(nnode,1);

    % Perform 4-point Gauss extrapolation to nodes
    xi_nodes = [-1,  1,  1, -1];
    eta_nodes = [-1, -1,  1,  1];

    for e = 1:mesh.nelem
        conn = mesh.conn(e,:);
        Xe   = mesh.nodes(conn,:);
        dofs = reshape([2*conn-1; 2*conn], [], 1);
        ue   = u(dofs);

        % Evaluate stress at 4 Gauss points
        gp_stress = struct();
        gp_stress.srr = zeros(4,1);
        gp_stress.stt = zeros(4,1);
        gp_stress.szz = zeros(4,1);
        gp_stress.srz = zeros(4,1);
        gp_stress.svm = zeros(4,1);

        for g = 1:mesh.ngp
            xi  = mesh.gp(g,1);
            eta = mesh.gp(g,2);
            [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
            [~, dNdX, ~]  = jacobian_2d(Xe, dNdxi);
            sigma = cauchy_stress_at_qp(Xe, ue, N, dNdX, par);

            srr = sigma(1,1);
            stt = sigma(2,2);
            szz = sigma(3,3);
            srz = sigma(1,3);

            gp_stress.srr(g) = srr;
            gp_stress.stt(g) = stt;
            gp_stress.szz(g) = szz;
            gp_stress.srz(g) = srz;
            gp_stress.svm(g) = sqrt(0.5*((srr-stt)^2 + (stt-szz)^2 + (szz-srr)^2 + 6*srz^2));
        end

        % Bilinear extrapolation matrix from Gauss points to nodes
        a = sqrt(3);
        for a_idx = 1:4
            xi_n = xi_nodes(a_idx);
            eta_n = eta_nodes(a_idx);

            E = 0.25 * [ (1 - a*xi_n)*(1 - a*eta_n), ...
                         (1 + a*xi_n)*(1 - a*eta_n), ...
                         (1 + a*xi_n)*(1 + a*eta_n), ...
                         (1 - a*xi_n)*(1 + a*eta_n) ];

            node = conn(a_idx);
            sigma_rr_sum(node) = sigma_rr_sum(node) + E * gp_stress.srr;
            sigma_tt_sum(node) = sigma_tt_sum(node) + E * gp_stress.stt;
            sigma_zz_sum(node) = sigma_zz_sum(node) + E * gp_stress.szz;
            sigma_rz_sum(node) = sigma_rz_sum(node) + E * gp_stress.srz;
            vm_sum(node)       = vm_sum(node)       + E * gp_stress.svm;
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


function plot_nodal_stress_contour(mesh, u, nodalStress, cbarLabel, plotTitle)
    nodes = mesh.nodes;
    conn  = mesh.conn;

    ur = u(1:2:end);
    uz = u(2:2:end);

    rDef = nodes(:,1) + ur;
    zDef = nodes(:,2) + uz;
    
    figure;

    patch('Faces', conn, ...
          'Vertices', [rDef*1e6, zDef*1e6], ...
          'FaceVertexCData', nodalStress, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');

    cb = colorbar;
    ylabel(cb, cbarLabel, 'Interpreter', 'tex');

    xlabel('r [\mum]');
    ylabel('z [\mum]');
    title(plotTitle, 'Interpreter', 'tex');
    xlim([0 15]);
    ylim([-1 5]);

    box on;
    set(gca,'FontSize',14);
end


function mesh = build_rounded_leukocyte_mesh(par)
    zvec = linspace(par.zMin, par.zMax, par.NzSolid).';

    rInner = par.RLin * ones(size(zvec));
    rOuter = rounded_outer_profile(zvec, par.RLout, par.zMin, par.Lz, par.Rc);

    nodes = zeros(par.NzSolid * par.NrL, 2);

    s = linspace(0,1,par.NrL);
    beta = 2;
    sBias = 1 - (1-s).^beta;

    for j = 1:par.NzSolid
        rline = rInner(j) + (rOuter(j)-rInner(j))*sBias;

        for i = 1:par.NrL
            id = sub2ind([par.NzSolid, par.NrL], j, i);
            nodes(id,:) = [rline(i), zvec(j)];
        end
    end

    conn = zeros((par.NrL-1)*(par.NzSolid-1),4);
    e = 0;

    for j = 1:par.NzSolid-1
        for i = 1:par.NrL-1
            n1 = sub2ind([par.NzSolid, par.NrL], j,   i);
            n2 = sub2ind([par.NzSolid, par.NrL], j,   i+1);
            n3 = sub2ind([par.NzSolid, par.NrL], j+1, i+1);
            n4 = sub2ind([par.NzSolid, par.NrL], j+1, i);
            e = e + 1;
            conn(e,:) = [n1 n2 n3 n4];
        end
    end

    gp1 = [-1, 1]/sqrt(3);
    gw1 = [1, 1];
    [g1,g2] = meshgrid(gp1,gp1);
    [w1,w2] = meshgrid(gw1,gw1);

    mesh = struct();
    mesh.nodes = nodes;
    mesh.conn  = conn;
    mesh.nelem = size(conn,1);
    mesh.ngp   = numel(g1);
    mesh.gp    = [g1(:), g2(:)];
    mesh.gw    = w1(:).*w2(:);
end


function plot_deformed_axisym_mesh(mesh, u, scale)
    if nargin < 3
        scale = 1;
    end

    r = mesh.nodes(:,1);
    z = mesh.nodes(:,2);

    ur = u(1:2:end);
    uz = u(2:2:end);

    rDef = r + scale*ur;
    zDef = z + scale*uz;

    patch('Faces', mesh.conn, ...
          'Vertices', [rDef*1e6, zDef*1e6], ...
          'FaceColor','none','EdgeColor','k');

    axis equal;
    xlabel('r [\mum]');
    ylabel('z [\mum]');
    box on;
end


function rOuter = rounded_outer_profile(z, Rout, zMin, L, Rc)
    zloc = z - zMin;

    rOuter = Rout * ones(size(zloc));

    idxL = zloc < Rc;
    xiL  = Rc - zloc(idxL);
    cutL = Rc - sqrt(max(Rc^2 - xiL.^2, 0));

    idxR = zloc > (L - Rc);
    xiR  = zloc(idxR) - (L - Rc);
    cutR = Rc - sqrt(max(Rc^2 - xiR.^2, 0));

    rOuter(idxL) = Rout - cutL;
    rOuter(idxR) = Rout - cutR;
end


