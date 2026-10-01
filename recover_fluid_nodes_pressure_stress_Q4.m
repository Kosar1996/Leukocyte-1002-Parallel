%function [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(meshF, ur2D, uz2D, mu,plist)
% before:
%function [pCell, sigmaCell, center] = recover_fluid_nodes_pressure_stress_Q4(meshF, ur2D, uz2D, mu, plist)

% after:
function [pCell, sigmaCell, center, pRaw] = recover_fluid_nodes_pressure_stress_Q4(meshF, ur2D, uz2D, mu, plist)
nodes = meshF.nodes;
elems = meshF.elems;
ne = size(elems,1);

pCell = zeros(ne,1);
%added_temp
pRaw = zeros(ne,1);
sigmaCell = zeros(ne,4); % [srr, stt, szz, srz]
center = zeros(ne,2);

for e = 1:ne
    conn = elems(e,:);
    xe = nodes(conn,1);
    ze = nodes(conn,2);

    [N, dNdxi] = shape_Q4(0,0);
    J = [xe ze].' * dNdxi;
    dNdx = dNdxi / J;
    dNdr = dNdx(:,1);
    dNdz = dNdx(:,2);

    r = N.' * xe;
    zc = N.' * ze;
    pc = N.' * plist(conn);

    ur = zeros(4,1);
    uz = zeros(4,1);
    for a = 1:4
        ur(a) = ur2D(conn(a));
        uz(a) = uz2D(conn(a));
    end

    durdr = dNdr.' * ur;
    duzdz = dNdz.' * uz;
    durdz = dNdz.' * ur;
    duzdr = dNdr.' * uz;
    ur_over_r = (N.' * ur) / r;

    divu = durdr + ur_over_r + duzdz;

    srr = -pc + 2*mu*durdr;
    stt = -pc + 2*mu*ur_over_r;
    szz = -pc + 2*mu*duzdz;
    srz = mu*(durdz + duzdr);

    pCell(e) = -(srr+stt+szz)/3;
    %added_temp
    pRaw(e) = pc;
    sigmaCell(e,:) = [srr, stt, szz, srz];
    center(e,:) = [r, zc];
end
end