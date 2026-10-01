function mesh = add_fluid_nodes(mesh)
nodeId = zeros(mesh.Nr,mesh.Nz);
nodes = zeros(mesh.Nr*mesh.Nz,2);

id = 0;

for j = 1:mesh.Nz
    for i = 1:mesh.Nr
        id = id + 1;
        nodeId(i,j) = id;
        nodes(id,:) = [mesh.Rp(i,j), mesh.Zp(i,j)];
    end
end

elems = zeros((mesh.Nr-1)*(mesh.Nz-1),4);
e = 0;
for j = 1:mesh.Nz-1
    for i = 1:mesh.Nr-1
        e = e + 1;
        n1 = nodeId(i,j);
        n2 = nodeId(i+1,j);
        n3 = nodeId(i+1,j+1);
        n4 = nodeId(i,j+1);
        elems(e,:) = [n1 n2 n3 n4];
    end
end

mesh.nodes = nodes;
mesh.elems = elems;
mesh.nodeId = nodeId;
end

