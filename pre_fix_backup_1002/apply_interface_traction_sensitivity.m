function [F, Kext, Bnormal, Btangent] = apply_interface_traction_sensitivity(mesh, u, F, interfaceNodes, traction)
% Robust traction sensitivity with respect to traction magnitudes on the
% fluid grid. This version is intentionally conservative: it reuses the
% existing apply_interface_traction(...) routine for the actual external
% force and follower-load tangent, then builds the force maps
% Bnormal = dFext/d(traction.normal) and Btangent = dFext/d(traction.tangent)
% by applying unit traction basis vectors.
%
% Why this replacement is needed:
% The previous implementation assumed that the number of traction samples
% equals the number of solid interface nodes. In the two-solid monolithic
% problem, traction.normal/tangent live on the fluid grid z, while the
% leukocyte mesh may have a different number of interface nodes. That caused
% "Index exceeds array bounds" when k+1 exceeded numel(traction.normal).
%
% This implementation works whether traction is defined on the reference
% fluid grid with traction.z or through traction_to_z(...).

    ndof = size(mesh.nodes,1)*2;
    Ntr = numel(traction.normal);

    % Actual force and follower-load tangent at the current traction.
    [F, Kext] = apply_interface_traction(mesh, u, F, interfaceNodes, traction);

    Bnormal  = sparse(ndof, Ntr);
    Btangent = sparse(ndof, Ntr);

    zeroTraction = traction;
    zeroTraction.normal  = zeros(size(traction.normal));
    zeroTraction.tangent = zeros(size(traction.tangent));

    for j = 1:Ntr
        trj = zeroTraction;
        trj.normal(j) = 1.0;
        Fj = zeros(ndof,1);
        Fj = apply_interface_traction(mesh, u, Fj, interfaceNodes, trj);
        Bnormal(:,j) = Fj;

        trj = zeroTraction;
        trj.tangent(j) = 1.0;
        Fj = zeros(ndof,1);
        Fj = apply_interface_traction(mesh, u, Fj, interfaceNodes, trj);
        Btangent(:,j) = Fj;
    end
end