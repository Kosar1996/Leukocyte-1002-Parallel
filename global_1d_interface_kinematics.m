function [rAtZ, UwAtZ, Hr, HU] = global_1d_interface_kinematics( ...
    mesh, uNew, uOld, interfaceNodes, zq, par)
%GLOBAL_1D_INTERFACE_KINEMATICS Embed finite solid interfaces in a larger 1D domain.
% Outside the actual solid-interface z range, the leukocyte side becomes the
% axis regularization radius and the endothelium side becomes the outer box
% radius. Those exterior rows have zero sensitivity to the solid DOFs.

    zq = zq(:);
    needSens = nargout > 2;
    ndof = 2 * size(mesh.nodes,1);
    isLeukocyte = min(mesh.nodes(:,1)) <= max(1e-12, 0.01 * max(par.RLout, 1e-12));
    if isLeukocyte
        exteriorRadius = global_1d_axis_radius(par);
    else
        exteriorRadius = global_1d_outer_radius(par);
    end

    ids = interfaceNodes(:);
    zNew = mesh.nodes(ids,2) + uNew(2*ids);
    zLo = min(zNew);
    zHi = max(zNew);
    tol = max(100 * eps(max(abs([zLo, zHi, zq(:).']))), 1e-12);
    active = find(zq >= zLo - tol & zq <= zHi + tol);
        
    rAtZ = exteriorRadius * ones(size(zq));
    UwAtZ = zeros(size(zq));

    if needSens
        Hr = sparse(numel(zq), ndof);
        HU = sparse(numel(zq), ndof);
    end

    if any(active)
        if needSens
            [rLoc, UwLoc, HrLoc, HULoc] = exact_interface_radius_velocity_sensitivity( ...
                mesh, uNew, uOld, interfaceNodes, zq(active), par.dt);
            Hr(active,:) = HrLoc;
            HU(active,:) = HULoc;
        else
            [rLoc, UwLoc] = exact_interface_radius_velocity_sensitivity( ...
                mesh, uNew, uOld, interfaceNodes, zq(active), par.dt);
        end
        rAtZ(active) = rLoc;
        UwAtZ(active) = UwLoc;
    end

    if isLeukocyte
        rAtZ = max(rAtZ, global_1d_axis_radius(par));
    else
        rAtZ = min(rAtZ, global_1d_outer_radius(par));
    end
end