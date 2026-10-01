function [rAtZ, UwAtZ, HrExact, HUExact, UrAtZ] = ...
    exact_interface_radius_velocity_sensitivity(mesh, uNew, uOld, interfaceNodes, zq, dt)
%EXACT_INTERFACE_RADIUS_VELOCITY_SENSITIVITY
% Exact deformed-interface mapping and first-order sensitivity.
%
% This evaluates the interface at fixed fluid-grid locations zq, while the
% solid interface itself is parameterized by its deformed axial coordinate
% z_def = Z + u_z. The derivative therefore includes the missing term
% dr(zq)/du_z caused by axial sliding of the deformed interface.
%
% Assumption: the deformed interface remains monotone in z, so the same
% bracketing segment is valid locally. If the sorted order changes during a
% finite-difference perturbation, the mapping is only piecewise smooth.
%
% UwAtZ is the AXIAL wall velocity only (d(u_z)/dt at fixed zq), matching
% the axial no-slip/slip boundary condition this model's governing
% equations actually use (lubrication-theory reduction: only u_z(r)
% enters the fluid BC, e.g. par draft eq. 8/10 "u_z - 1 = -l*du_z/dr").
% The interface's radial motion enters the fluid problem separately,
% through the evolving boundary radius rAtZ itself (via the continuity/
% evolution equation, draft eq. 11) -- NOT through a wall-velocity BC.
% For a full 3D velocity, UrAtZ (radial
% interface velocity, d(r)/dt at fixed zq) is now also returned as an
% additional, backward-compatible output (existing callers requesting
% <=4 outputs are unaffected). It is not currently consumed by any
% boundary condition -- see the accompanying report for why, and flag if
% a specific new BC should use it.

    ids = interfaceNodes(:);
    ndof = 2 * size(mesh.nodes,1);

    rNew  = mesh.nodes(ids,1) + uNew(2*ids - 1);
    zNew  = mesh.nodes(ids,2) + uNew(2*ids);
    uzNew = uNew(2*ids);

    rOld  = mesh.nodes(ids,1) + uOld(2*ids - 1);
    zOld  = mesh.nodes(ids,2) + uOld(2*ids);
    uzOld = uOld(2*ids);

    [zNewS, idxNew] = sort(zNew);
    rNewS  = rNew(idxNew);
    uzNewS = uzNew(idxNew);
    idsNew = ids(idxNew);

    [zOldS, idxOld] = sort(zOld);
    uzOldS = uzOld(idxOld);
    rOldS  = rOld(idxOld);

    zq = zq(:);
    nq = numel(zq);
    ns = numel(zNewS);

    if ns < 2
        error('Need at least two interface nodes for exact-interface interpolation.');
    end

    rAtZ = zeros(nq,1);
    uzNewAtZ = zeros(nq,1);

    if nargout > 2
        rows = zeros(4*nq,1);
        cols = zeros(4*nq,1);
        valsR = zeros(4*nq,1);
        valsU = zeros(4*nq,1);
        ptr = 1;
    end

    for q = 1:nq
        zc = zq(q);

        if zc <= zNewS(1)
            k = 1;
        elseif zc >= zNewS(end)
            k = ns - 1;
        else
            k = find(zNewS <= zc, 1, 'last');
            if k >= ns
                k = ns - 1;
            end
        end

        z1 = zNewS(k);
        z2 = zNewS(k+1);
        r1 = rNewS(k);
        r2 = rNewS(k+1);
        w1 = uzNewS(k);
        w2 = uzNewS(k+1);

        id1 = idsNew(k);
        id2 = idsNew(k+1);

        L = z2 - z1;
        if abs(L) < 100*eps(max(1, max(abs(zNewS))))
            error('Degenerate deformed interface segment in z.');
        end

        a = (zc - z1) / L;

        rAtZ(q)     = (1-a)*r1 + a*r2;
        uzNewAtZ(q) = (1-a)*w1 + a*w2;

        if nargout > 2
            % a = (zq - z1)/(z2-z1)
            da_dz1 = (zc - z2) / L^2;
            da_dz2 = -(zc - z1) / L^2;

            % r(zq) = (1-a) r1 + a r2
            dr_dur1 = 1 - a;
            dr_dur2 = a;
            dr_duz1 = (r2 - r1) * da_dz1;
            dr_duz2 = (r2 - r1) * da_dz2;

            % w(zq) = (1-a) w1 + a w2, with z1,z2 depending on w1,w2.
            dw_duz1 = (1 - a) + (w2 - w1) * da_dz1;
            dw_duz2 = a       + (w2 - w1) * da_dz2;

            nodeTol = 100*eps(max(1, max(abs(zNewS))));
            if abs(zc - z1) <= nodeTol && k > 1
                Lleft = zNewS(k) - zNewS(k-1);
                if abs(Lleft) >= 100*eps(max(1, max(abs(zNewS))))
                    rSlopeLeft = (rNewS(k) - rNewS(k-1)) / Lleft;
                    wSlopeLeft = (uzNewS(k) - uzNewS(k-1)) / Lleft;
                    rSlopeRight = (r2 - r1) / L;
                    wSlopeRight = (w2 - w1) / L;
                    dr_duz1 = -0.5 * (rSlopeLeft + rSlopeRight);
                    dw_duz1 = 1 - 0.5 * (wSlopeLeft + wSlopeRight);
                end
            elseif abs(zc - z2) <= nodeTol && k < ns - 1
                Lright = zNewS(k+2) - zNewS(k+1);
                if abs(Lright) >= 100*eps(max(1, max(abs(zNewS))))
                    rSlopeLeft = (r2 - r1) / L;
                    wSlopeLeft = (w2 - w1) / L;
                    rSlopeRight = (rNewS(k+2) - rNewS(k+1)) / Lright;
                    wSlopeRight = (uzNewS(k+2) - uzNewS(k+1)) / Lright;
                    dr_duz2 = -0.5 * (rSlopeLeft + rSlopeRight);
                    dw_duz2 = 1 - 0.5 * (wSlopeLeft + wSlopeRight);
                end
            end

            rows(ptr:ptr+3) = q;
            cols(ptr:ptr+3) = [2*id1-1; 2*id2-1; 2*id1; 2*id2];
            valsR(ptr:ptr+3) = [dr_dur1; dr_dur2; dr_duz1; dr_duz2];
            valsU(ptr:ptr+3) = [0; 0; dw_duz1/dt; dw_duz2/dt];
            ptr = ptr + 4;
        end
    end

    uzOldAtZ = interp_curve_values(zOldS, uzOldS, zq);
    UwAtZ = (uzNewAtZ - uzOldAtZ) / dt;

    if nargout > 4
        % Radial interface velocity at the same fixed zq query points,
        % same construction as UwAtZ above: interpolate the OLD state's
        % radius at zq (using the old deformed curve), subtract from the
        % already-computed NEW radius (rAtZ), divide by dt.
        rOldAtZ = interp_curve_values(zOldS, rOldS, zq);
        UrAtZ = (rAtZ - rOldAtZ) / dt;
    end

    if nargout > 2
        HrExact = sparse(rows, cols, valsR, nq, ndof);
        HUExact = sparse(rows, cols, valsU, nq, ndof);
    end
end