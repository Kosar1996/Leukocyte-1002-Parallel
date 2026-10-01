function [rAtZ, wzAtZ, wrAtZ] = exact_interface_radius_velocity( ...
    mesh, uNew, uOld, interfaceNodes, zq, dt,par)
% EXACT_INTERFACE_RADIUS_VELOCITY
% Position-aware interpolation of deformed interface radius and flux-equivalent 
% radial/axial mesh boundary velocities with domain-aware boundary clamping.
    isLeukocyte = min(mesh.nodes(:,1)) <= max(1e-12, 0.01 * max(par.RLout, 1e-12));
    if isLeukocyte
        exteriorRadius = global_1d_axis_radius(par);
    else
        exteriorRadius = global_1d_outer_radius(par);
    end
    % Extract deformed coordinates and displacements
    [rNew, zNew, uzNew] = deformed_interface_curve(mesh, uNew, interfaceNodes);
    [rOld,    zOld, uzOld] = deformed_interface_curve(mesh, uOld, interfaceNodes);

    % Interpolate current and previous radius along current spatial curve
    rAtZ    = interp_curve_values(zNew, rNew, zq);
    rOldAtZ = interp_curve_values(zOld, rOld, zq);

    % % Identify active axial bounds of the mesh interface
    % zMinNew = min(zNew);
    % zMaxNew = max(zNew);
        % Endpoint clamping safeguard to prevent artificial boundary bulges near zMin and zMax
    rAtZ(zq < min(zNew)) = exteriorRadius;
    rAtZ(zq > max(zNew)) = exteriorRadius;
        rOldAtZ(zq < min(zOld)) = exteriorRadius;
    rOldAtZ(zq > max(zOld)) = exteriorRadius;
     isOut   = (zq < min(zNew) )| (zq > max(zNew));


    % Interpolate axial displacement field
    uzNewAtZ = interp_curve_values(zNew, uzNew, zq);
    uzOldAtZ = interp_curve_values(zOld, uzOld, zq);

    % Zero axial motion outside the physical solid body extent
    uzNewAtZ(isOut) = 0;
    isOut   = (zq < min(zOld) )| (zq > max(zOld));
    uzOldAtZ(isOut) = 0;

    % Material axial mesh interface velocity
    wzAtZ = (uzNewAtZ - uzOldAtZ) / dt;

% 1. Geometric quadratic radial sweep rate
    rMidAtZ = 0.5 * (rAtZ + rOldAtZ);
    wr_geometric = (rAtZ.^2 - rOldAtZ.^2) ./ (2 * max(rMidAtZ, eps) * dt);

    % 2. Surface slope term (dr/dz) along the target grid zq
    dzq = gradient(zq(:));
    dr_dz = gradient(rAtZ(:)) ./ max(dzq, eps);

    % 3. Effective radial wall velocity incorporating axial mesh motion along slope
    wrAtZ = wr_geometric(:) - wzAtZ(:) .* dr_dz;

    % Zero velocities outside physical solid extent
         isOutNew   = (zq < min(zNew) )| (zq > max(zNew));
    wrAtZ(isOutNew) = 0;
    wzAtZ(isOutNew) = 0;
end