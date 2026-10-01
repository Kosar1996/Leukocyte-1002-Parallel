function Fint = assemble_finite_def_internal_force_only(mesh, u, par)

    ndof = size(mesh.nodes,1)*2;
    Fint = zeros(ndof,1);
    useCache = isfield(mesh, 'axisymCache');
    if useCache
        cache = mesh.axisymCache;
    end

    for e = 1:mesh.nelem
        if useCache
            dofs = cache.dofs(e,:).';
            fe = finite_def_element_residual_only_cached(cache, e, u(dofs), par);
        else
            conn = mesh.conn(e,:);
            Xe   = mesh.nodes(conn,:);
            dofs = reshape([2*conn-1; 2*conn], [], 1);
            fe = finite_def_element_residual_only(Xe, u(dofs), mesh, par);
        end
        Fint(dofs) = Fint(dofs) + fe;
    end
end

function fe = finite_def_element_residual_only_cached(cache, e, ue, par)
    fe = zeros(8,1);

    Rnod = cache.Rnod(:,e);
    Znod = cache.Znod(:,e);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);
    I3 = eye(3);

    for g = 1:cache.ngp
        N = cache.N(:,g,e).';
        dNdX = cache.dNdX(:,:,g,e);
        detJ0 = cache.detJ0(g,e);
        Rg = cache.Rg0(g,e);
        w = cache.gw(g);

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
        trB = B(1,1) + B(2,2) + B(3,3);
        T = par.Ge * J^(-5/3) * ( B - (trB/3)*I3 ) ...
          + par.Ke * (J - 1) * I3;

        P = (F \ (J * T).').';
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

function fe = finite_def_element_residual_only(Xe, ue, mesh, par)

    fe = zeros(8,1);

    Rnod = Xe(:,1);
    Znod = Xe(:,2);

    rnod = Rnod + ue(1:2:end);
    znod = Znod + ue(2:2:end);

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
        I = eye(3);
        trB = B(1,1) + B(2,2) + B(3,3);

        % Exact constitutive law requested by user
        T = par.Ge * J^(-5/3) * ( B - (trB/3)*I ) ...
          + par.Ke * (J - 1) * I;

        P = (F \ (J * T).').';

        for a = 1:4
            dNa_dR = dNdX(a,1);
            dNa_dZ = dNdX(a,2);
            Na     = N(a);

            fe(2*a-1) = fe(2*a-1) + ...
                ( P(1,1)*dNa_dR + P(1,3)*dNa_dZ + P(2,2)*(Na/Rg) ) ...
                * (2*pi*Rg) * detJ0 * w;

            fe(2*a) = fe(2*a) + ...
                ( P(3,1)*dNa_dR + P(3,3)*dNa_dZ ) ...
                * (2*pi*Rg) * detJ0 * w;
        end
    end
end