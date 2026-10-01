function [rDef, zDef] = local_deformed_interface_curve(mesh, u, interfaceNodes)
    ids = interfaceNodes(:);
    rDef = mesh.nodes(ids,1) + u(2*ids - 1);
    zDef = mesh.nodes(ids,2) + u(2*ids);
    [zDef, idx] = sort(zDef);
    rDef = rDef(idx);
end
