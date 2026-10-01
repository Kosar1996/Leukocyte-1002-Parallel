function solidMask = deformed_solid_mask_on_grid(mesh, u, r, z)
    solidMask = false(numel(r), numel(z));
    if isempty(mesh) || isempty(u) || ~isfield(mesh, 'nodes') || ...
            ~isfield(mesh, 'conn') || numel(u) < 2*size(mesh.nodes,1)
        return;
    end

    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    [Rgrid, Zgrid] = ndgrid(r(:), z(:));

    for e = 1:size(mesh.conn,1)
        ids = mesh.conn(e,:);
        rv = rDef(ids);
        zv = zDef(ids);
        inBox = Rgrid >= min(rv) & Rgrid <= max(rv) & ...
                Zgrid >= min(zv) & Zgrid <= max(zv);
        if any(inBox(:))
            inLocal = false(size(solidMask));
            inLocal(inBox) = inpolygon(Rgrid(inBox), Zgrid(inBox), rv, zv);
            solidMask = solidMask | inLocal;
        end
    end
end