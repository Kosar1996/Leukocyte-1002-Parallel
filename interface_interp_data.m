function data = interface_interp_data(mesh, interfaceNodes, zq)
    persistent cache

    interfaceNodes = interfaceNodes(:);
    zq = zq(:);
    nnode = size(mesh.nodes,1);

    if isempty(cache)
        cache = struct('nnode', {}, 'interfaceNodes', {}, 'zq', {}, ...
            'Hr', {}, 'Hz', {}, 'rBase', {});
    end

    for kCache = 1:numel(cache)
        if cache(kCache).nnode == nnode && ...
                isequal(cache(kCache).interfaceNodes, interfaceNodes) && ...
                isequal(cache(kCache).zq, zq)
            data = cache(kCache);
            return;
        end
    end

    ndof = nnode * 2;
    Nz = numel(zq);
    zn = mesh.nodes(interfaceNodes,2);
    [zs, idx] = sort(zn);
    ids = interfaceNodes(idx);

    rows = zeros(2*Nz,1);
    colsR = zeros(2*Nz,1);
    colsZ = zeros(2*Nz,1);
    vals = zeros(2*Nz,1);
    rBase = zeros(Nz,1);

    for j = 1:Nz
        zj = zq(j);

        if zj <= zs(1)
            iSeg = 1;
            w2 = 0;
        elseif zj >= zs(end)
            iSeg = numel(zs) - 1;
            w2 = 1;
        else
            iSeg = find(zs <= zj, 1, 'last');
            if iSeg == numel(zs)
                iSeg = iSeg - 1;
            end
            w2 = (zj - zs(iSeg)) / (zs(iSeg+1) - zs(iSeg));
        end

        w1 = 1 - w2;
        n1 = ids(iSeg);
        n2 = ids(iSeg+1);
        loc = (2*j-1):(2*j);

        rows(loc) = j;
        colsR(loc) = [2*n1-1; 2*n2-1];
        colsZ(loc) = [2*n1; 2*n2];
        vals(loc) = [w1; w2];
        rBase(j) = w1 * mesh.nodes(n1,1) + w2 * mesh.nodes(n2,1);
    end

    data = struct();
    data.nnode = nnode;
    data.interfaceNodes = interfaceNodes;
    data.zq = zq;
    data.Hr = sparse(rows, colsR, vals, Nz, ndof);
    data.Hz = sparse(rows, colsZ, vals, Nz, ndof);
    data.rBase = rBase;

    cache(end+1) = data;
    if numel(cache) > 8
        cache = cache(end-7:end);
    end
end