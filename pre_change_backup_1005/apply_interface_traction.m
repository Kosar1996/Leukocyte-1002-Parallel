function [F, Kext] = apply_interface_traction(mesh, u, F, interfaceNodes, traction)
% APPLY_INTERFACE_TRACTION
% Adds distributed follower traction on the deformed interface to the global 
% external force vector F, and returns the consistent follower tangent Kext.

    ndof = size(mesh.nodes,1) * 2;
    Kext = sparse(ndof, ndof);

    % Preserve topological node ordering directly from mesh definition
    interface = interfaceNodes(:);
    
    % Evaluate reference and deformed axial coordinates along ordered nodes
    zn_ref = mesh.nodes(interface, 2);
    zn_def = zn_ref + u(2*interface);

    % Interpolate traction magnitudes along appropriate spatial configuration
    if isfield(traction, 'z')
        tr_n = interp_curve_values(traction.z, traction.normal, zn_def);
        tr_t = interp_curve_values(traction.z, traction.tangent, zn_def);
        [tr_n, tr_t] = apply_traction_support_window(traction, zn_def, tr_n, tr_t);
    else
        zq = traction_to_z(traction, zn_ref);
        tr_n = zq.normal;
        tr_t = zq.tangent;
        [tr_n, tr_t] = apply_traction_support_window(traction, zn_ref, tr_n, tr_t);
    end

    % -------------------------------------------------------------------------
    % CORRECTED ORIENTATION SIGN LOGIC:
    % Outer boundary (Endothelium): normalSign = +1.0 -> pushes +r into tissue
    % Inner boundary (Leukocyte)  : normalSign = -1.0 -> pushes -r into cell core
    % -------------------------------------------------------------------------
    isInnerBoundary = isfield(traction, 'isInnerBoundary') && traction.isInnerBoundary;
    normalSign = 1.0; 
    if isInnerBoundary
        normalSign = -1.0;
    end

    % 2-point Gauss quadrature rule
    xi_gp = [-1, 1] / sqrt(3);
    w_gp  = [1, 1];

    for k = 1:numel(interface)-1
        n1 = interface(k);
        n2 = interface(k+1);

        dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

        % Reference coordinates
        r1 = mesh.nodes(n1,1);  z1 = mesh.nodes(n1,2);
        r2 = mesh.nodes(n2,1);  z2 = mesh.nodes(n2,2);

        % Current deformed coordinates
        ur1 = u(2*n1-1); uz1 = u(2*n1);
        ur2 = u(2*n2-1); uz2 = u(2*n2);

        x1 = [r1 + ur1; z1 + uz1];
        x2 = [r2 + ur2; z2 + uz2];

        dx = x2 - x1;
        L  = norm(dx);

        if L <= 1e-14
            continue;
        end

        t_hat = dx / L;

        % Unit normal vector with corrected boundary orientation
        n_hat = normalSign * [ t_hat(2); -t_hat(1) ];

        % Nodal traction magnitudes
        tn_nodes = [tr_n(k);   tr_n(k+1)];
        tt_nodes = [tr_t(k);   tr_t(k+1)];

        fe = zeros(4,1);
        ke = zeros(4,4);

        for g = 1:2
            xi = xi_gp(g);
            wg = w_gp(g);

            N1 = 0.5 * (1 - xi);
            N2 = 0.5 * (1 + xi);

            Nline = [N1, N2];

            % Scalar radius at Gauss point in current configuration
            r_gp = N1 * x1(1) + N2 * x2(1);

            tn_gp = Nline * tn_nodes;
            tt_gp = Nline * tt_nodes;

            % Current traction vector in spatial basis
            tvec = tn_gp * n_hat + tt_gp * t_hat;

            Nmat = [N1 0  N2 0;
                    0  N1 0  N2];

            Jline = L / 2;
            fac   = (2*pi*r_gp) * Jline * wg;

            % External force vector contribution
            fe = fe + (Nmat.' * tvec) * fac;

            % Consistent follower load derivatives
            Bdx = [-1  0  1  0;
                    0 -1  0  1];

            Br  = [N1 0 N2 0];

            I2 = eye(2);
            Ptan = I2 - (t_hat * t_hat.');

            % Rotation operator matching boundary normal orientation
            R90 = normalSign * [0 1; -1 0];

            At = (Ptan / L) * Bdx;
            An = R90 * At;

            AJ = 0.5 * (t_hat.' * Bdx);
            Afac = 2*pi * wg * ( Jline * Br + r_gp * AJ );
            Avec = tn_gp * An + tt_gp * At;

            % Consistent element tangent matrix
            ke = ke + (Nmat.' * Avec) * fac + (Nmat.' * tvec) * Afac;
        end

        F(dofs) = F(dofs) + fe;
        Kext(dofs,dofs) = Kext(dofs,dofs) + ke;
    end
end

function zq = traction_to_z(traction, zNodes)
    N = numel(traction.normal);
    z0 = linspace(min(zNodes), max(zNodes), N).';
    zq.normal  = safe_interp1_same_or_resample(z0, traction.normal,  zNodes, 'traction.normal');
    zq.tangent = safe_interp1_same_or_resample(z0, traction.tangent, zNodes, 'traction.tangent');
end

function [tr_n, tr_t] = apply_traction_support_window(traction, zs, tr_n, tr_t)
    useWindow = isfield(traction, 'outsideZero') && traction.outsideZero;
    if ~useWindow
        return;
    end

    if isfield(traction, 'supportInterval') && numel(traction.supportInterval) == 2
        support = sort(traction.supportInterval(:));
    elseif isfield(traction, 'z') && ~isempty(traction.z)
        support = [min(traction.z(:)); max(traction.z(:))];
    elseif isempty(zs)
        support = [0; 0];
    else
        support = [min(zs(:)); max(zs(:))];
    end

    outside = zs(:) < support(1) | zs(:) > support(2);
    tr_n(outside) = 0;
    tr_t(outside) = 0;
end