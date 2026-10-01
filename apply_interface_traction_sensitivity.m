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

    % SPEED-UP (Oct 2, results bit-identical to the previous version):
    %  - the 2*Ntr unit-traction calls only need the force vector, so they are
    %    called with one output (no follower tangent is built, see
    %    apply_interface_traction.m);
    %  - columns are collected in dense arrays and converted to sparse once
    %    (same values/pattern as assigning Bnormal(:,j) = Fj column by column);
    %  - the force-only calls use a vectorized path (all segments at once, same
    %    arithmetic, see apply_interface_traction.m);
    %  - optional parfor over the 2*Ntr independent calls (TRACTION_PARFOR=1);
    %    off by default because the vectorized calls are faster than the
    %    parfor overhead. Only mesh.nodes is sent to the workers.
    zeroTraction = traction;
    zeroTraction.normal  = zeros(size(traction.normal));
    zeroTraction.tangent = zeros(size(traction.tangent));
    meshNodesOnly = struct('nodes', mesh.nodes);

    BnD = zeros(ndof, Ntr);
    BtD = zeros(ndof, Ntr);
    if use_traction_parfor()
        parfor j = 1:Ntr
            [BnD(:,j), BtD(:,j)] = unit_traction_columns(meshNodesOnly, u, ndof, interfaceNodes, zeroTraction, j);
        end
    else
        for j = 1:Ntr
            [BnD(:,j), BtD(:,j)] = unit_traction_columns(meshNodesOnly, u, ndof, interfaceNodes, zeroTraction, j);
        end
    end
    Bnormal  = sparse(BnD);
    Btangent = sparse(BtD);
end

function [fn, ft] = unit_traction_columns(mesh, u, ndof, interfaceNodes, zeroTraction, j)
    trj = zeroTraction;
    trj.normal(j) = 1.0;
    fn = zeros(ndof,1);
    fn = apply_interface_traction(mesh, u, fn, interfaceNodes, trj);

    trj = zeroTraction;
    trj.tangent(j) = 1.0;
    ft = zeros(ndof,1);
    ft = apply_interface_traction(mesh, u, ft, interfaceNodes, trj);
end

function tf = use_traction_parfor()
% parfor over the unit-traction columns is OFF by default: with the
% vectorized force-only path one full sensitivity call takes ~0.02 s, less
% than the parfor dispatch/transfer overhead. TRACTION_PARFOR=1 switches it
% on (only if a pool is already open; never opens a pool itself).
    persistent envOn
    if isempty(envOn), envOn = strcmp(getenv('TRACTION_PARFOR'), '1'); end
    tf = false;
    if ~envOn, return; end
    try
        tf = ~isempty(gcp('nocreate'));
    catch
        tf = false;
    end
end
