% CHANGES TO LOOK FOR IN THIS FILE:
% - Lines 12-24 (Issue 1): Added explicit dt fallback check to prevent Division-by-Zero 
%   (NaN/Inf) in wall velocity (Uw) calculations when par.dt is zero, unpopulated, or missing.
% - Lines 110-116 (Issue 1): Guarded UwCand evaluation in physical_endothelium_fluid_boundary_kinematics
%   to return zero velocities during dt == 0 diagnostic or warm-start passes.
% - Lines 120-134 (Issue 2): Sanitized envelope selection in physical_endothelium_fluid_boundary_kinematics 
%   by validating radius candidates (rCand > 0 & isfinite) before evaluating minimum boundary updates, 
%   preventing internal mesh edges or re-entrant cuts from corrupting the fluid-facing wall profile.

function state = attach_physical_solid_interface_fields( ...
    state, old, meshE, interfaceE, meshL, interfaceL, z, par)
% Store the actual FE solid interfaces separately from global reservoir
% embeddings such as r = 15 um outside the endothelium span.

    % OLD BUGGY CODE (Issue 1): Passed par.dt directly without checking if it was populated or > 0:
    % [state.deltaESolid, state.UwESolid, ...] = physical_solid_interface_kinematics(meshE, state.uE, old.uE, interfaceE, z, par.dt);

    % FIXED (Issue 1): Sanitize time step dt to prevent division-by-zero across sub-functions
    dt = 1.0;
    zeroUw = true;
    if isfield(par, 'dt') && isfinite(par.dt) && par.dt > 0
        dt = par.dt;
        zeroUw = false;
    end

    [state.deltaESolid, state.UwESolid, ...
        state.endotheliumInterfaceActive, state.endotheliumInterfaceZRange] = ...
        physical_solid_interface_kinematics( ...
            meshE, state.uE, old.uE, interfaceE, z, dt);
    [state.deltaEFluidBoundary, state.UwEFluidBoundary, ...
        state.endotheliumFluidBoundaryActive, state.endotheliumFluidBoundaryZRange] = ...
        physical_endothelium_fluid_boundary_kinematics( ...
            meshE, state.uE, old.uE, z, dt);

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        state.deltaLSolid = state.deltaL(:);
        state.UwLSolid = zeros(size(z(:)));
        state.leukocyteInterfaceActive = true(size(z(:)));
        state.leukocyteInterfaceZRange = [min(z(:)); max(z(:))];
    elseif use_RLout_fluid_interface_for_solid_leukocyte(par) || ...
            isempty(meshL) || isempty(interfaceL) || ~isfield(state, 'uL') || isempty(state.uL)
        state.deltaLSolid = state.deltaL(:);
        state.UwLSolid = state.UwL(:);
        state.leukocyteInterfaceActive = true(size(z(:)));
        state.leukocyteInterfaceZRange = [min(z(:)); max(z(:))];
    else
        [state.deltaLSolid, state.UwLSolid, ...
            state.leukocyteInterfaceActive, state.leukocyteInterfaceZRange] = ...
            physical_solid_interface_kinematics( ...
                meshL, state.uL, old.uL, interfaceL, z, dt);
    end

    % Zero out velocities if dt was unpopulated or zero
    if zeroUw
        state.UwESolid(:) = 0;
        state.UwEFluidBoundary(:) = 0;
        state.UwLSolid(:) = 0;
    end
end

function [rAtZ, UwAtZ, active, zRange] = physical_solid_interface_kinematics( ...
    mesh, uNew, uOld, interfaceNodes, zq, dt)
    zq = zq(:);
    [rAtZ, UwAtZ] = exact_interface_radius_velocity_sensitivity( ...
        mesh, uNew, uOld, interfaceNodes, zq, dt);

    ids = interfaceNodes(:);
    zDef = mesh.nodes(ids, 2) + uNew(2*ids);
    zRange = [min(zDef); max(zDef)];
    tol = max(100 * eps(max(abs([zRange(:); zq(:); 1]))), 1e-12);
    active = zq >= zRange(1) - tol & zq <= zRange(2) + tol;
end

function [rAtZ, UwAtZ, active, zRange] = physical_endothelium_fluid_boundary_kinematics( ...
    mesh, uNew, uOld, zq, dt)
% Envelope of the full fluid-facing endothelium boundary.
    zq = zq(:);
    edges = endothelium_fluid_boundary_edges(mesh);
    if isempty(edges)
        rAtZ = nan(size(zq));
        UwAtZ = zeros(size(zq));
        active = false(size(zq));
        zRange = [NaN; NaN];
        return;
    end

    rNew = mesh.nodes(:,1) + uNew(1:2:end);
    zNew = mesh.nodes(:,2) + uNew(2:2:end);
    uzNew = uNew(2:2:end);
    uzOld = uOld(2:2:end);

    zEdge = zNew(edges(:));
    zEdge = zEdge(isfinite(zEdge));
    if isempty(zEdge)
        rAtZ = nan(size(zq));
        UwAtZ = zeros(size(zq));
        active = false(size(zq));
        zRange = [NaN; NaN];
        return;
    end
    zRange = [min(zEdge); max(zEdge)];

    rAtZ = inf(size(zq));
    UwAtZ = zeros(size(zq));
    tol = max(100 * eps(max(abs([zRange(:); zq(:); 1]))), 1e-12);

    for e = 1:size(edges,1)
        n1 = edges(e,1);
        n2 = edges(e,2);
        z1 = zNew(n1);
        z2 = zNew(n2);
        r1 = rNew(n1);
        r2 = rNew(n2);
        if ~all(isfinite([z1, z2, r1, r2]))
            continue;
        end

        zLo = min(z1, z2);
        zHi = max(z1, z2);
        onSegment = zq >= zLo - tol & zq <= zHi + tol;
        if ~any(onSegment)
            continue;
        end

        if abs(z2 - z1) <= tol
            alpha = 0.5 * ones(nnz(onSegment), 1);
        else
            alpha = (zq(onSegment) - z1) / (z2 - z1);
            alpha = min(1, max(0, alpha));
        end
        rCand = (1 - alpha) * r1 + alpha * r2;
        uzNewCand = (1 - alpha) * uzNew(n1) + alpha * uzNew(n2);
        uzOldCand = (1 - alpha) * uzOld(n1) + alpha * uzOld(n2);

        % OLD BUGGY LINE (Issue 1): UwCand = (uzNewCand - uzOldCand) / dt;
        % FIXED (Issue 1): Explicit dt check prevents Inf/NaN propagation when dt == 0
        if dt > 0
            UwCand = (uzNewCand - uzOldCand) / dt;
        else
            UwCand = zeros(size(uzNewCand));
        end

        ids = find(onSegment);

        % OLD BUGGY CODE (Issue 2): Unconditionally evaluated rCand < rAtZ(ids), allowing NaN/negative
        % or internal mesh cutout edges to corrupt the outer fluid boundary envelope:
        % improve = rCand < rAtZ(ids);
        % rAtZ(ids(improve)) = rCand(improve);
        % UwAtZ(ids(improve)) = UwCand(improve);

        % FIXED (Issue 2): Explicitly validate radius candidates (rCand > 0 & isfinite) before 
        % updating the fluid boundary envelope:
        validCand = isfinite(rCand) & (rCand > 0);
        idsValid = ids(validCand);
        rValid = rCand(validCand);
        UwValid = UwCand(validCand);

        improve = rValid < rAtZ(idsValid);
        rAtZ(idsValid(improve)) = rValid(improve);
        UwAtZ(idsValid(improve)) = UwValid(improve);
    end

    active = isfinite(rAtZ);
    rAtZ(~active) = NaN;
end

function edges = endothelium_fluid_boundary_edges(mesh)
    edges = boundary_edges_from_q4_mesh(mesh.conn);
    if isempty(edges)
        return;
    end
    rRef = mesh.nodes(:,1);
    rOuter = max(rRef);
    tol = max(1e-12, 1e-8 * max(abs(rOuter), realmin));
    outerEdge = abs(rRef(edges(:,1)) - rOuter) <= tol & ...
        abs(rRef(edges(:,2)) - rOuter) <= tol;
    edges = edges(~outerEdge,:);
end

function edges = boundary_edges_from_q4_mesh(conn)
    allEdges = [conn(:,[1 2]); conn(:,[2 3]); conn(:,[3 4]); conn(:,[4 1])];
    allEdges = sort(allEdges, 2);
    [edgeUnique, ~, ic] = unique(allEdges, 'rows');
    counts = accumarray(ic, 1);
    edges = edgeUnique(counts == 1, :);
end