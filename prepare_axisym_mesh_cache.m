function mesh = prepare_axisym_mesh_cache(mesh)
    if isfield(mesh, 'axisymCache') && ...
            isfield(mesh.axisymCache, 'version') && mesh.axisymCache.version == 1
        return;
    end

    nelem = mesh.nelem;
    ngp = mesh.ngp;

    dofs = zeros(nelem, 8);
    XeAll = zeros(4, 2, nelem);
    RnodAll = zeros(4, nelem);
    ZnodAll = zeros(4, nelem);
    Nall = zeros(4, ngp, nelem);
    dNdXall = zeros(4, 2, ngp, nelem);
    detJ0all = zeros(ngp, nelem);
    Rg0all = zeros(ngp, nelem);
    iK = zeros(nelem * 64, 1);
    jK = zeros(nelem * 64, 1);

    for e = 1:nelem
        conn = mesh.conn(e,:);
        Xe = mesh.nodes(conn,:);
        dofRow = reshape([2*conn-1; 2*conn], [], 1).';

        dofs(e,:) = dofRow;
        XeAll(:,:,e) = Xe;
        RnodAll(:,e) = Xe(:,1);
        ZnodAll(:,e) = Xe(:,2);

        [ii, jj] = ndgrid(dofRow, dofRow);
        loc = (64*(e-1)+1):(64*e);
        iK(loc) = ii(:);
        jK(loc) = jj(:);

        for g = 1:ngp
            xi = mesh.gp(g,1);
            eta = mesh.gp(g,2);
            [N, dNdxi, ~] = q4_shape(xi, eta, 1.0);
            [~, dNdX, detJ0] = jacobian_2d(Xe, dNdxi);

            Nall(:,g,e) = N(:);
            dNdXall(:,:,g,e) = dNdX;
            detJ0all(g,e) = detJ0;
            Rg0all(g,e) = N * Xe(:,1);
        end
    end

    mesh.axisymCache = struct( ...
        'version', 1, ...
        'nelem', nelem, ...
        'ngp', ngp, ...
        'gw', mesh.gw(:), ...
        'dofs', dofs, ...
        'Xe', XeAll, ...
        'Rnod', RnodAll, ...
        'Znod', ZnodAll, ...
        'N', Nall, ...
        'dNdX', dNdXall, ...
        'detJ0', detJ0all, ...
        'Rg0', Rg0all, ...
        'iK', iK, ...
        'jK', jK);
end