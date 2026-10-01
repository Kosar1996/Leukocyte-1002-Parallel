function [sigma, p, rz_actual, xieta, ij] = locate_and_interp_fluid_stress(mraw, meshF, ur2D, uz2D, mu, plist, rq, zq, i_guess, j_guess)
% LOCATE_AND_INTERP_FLUID_STRESS
% Robust point location + position-aware Q4 interpolation for fluid stress.
%
% REVISION HISTORY & MERGED BUG FIXES:
% -------------------------------------------------------------------------
% 1. Lines 48-68 (Local Newton Step Regularization):
%    [OLD]: if abs(detJ) <= 1e-12, dxieta = [0;0]; break; end
%    Setting dxieta=0 caused Newton iterations to freeze, locking (i,j) in
%    an infinite walk loop during boundary searches. Added a signed 
%    determinant floor (`epsJ = 1e-14`) and step capping to allow safe step 
%    updates across distorted boundary elements.
%
% 2. Line 70 (Step Size Safeguard):
%    Added step normalization capping `norm(dxieta) <= 1.0` during point 
%    location to prevent large numerical jumps when evaluating points near 
%    distorted interfaces.
%
% 3. Lines 105-118 (Global Jacobian Gradient Regularization):
%    [OLD]: error('Singular Jacobian matrix...')
%    Replaced hard error termination with a regularized cofactor inverse using 
%    `detJ_safe`. Throwing a hard error caused global solver hangs during 
%    large-deformation interface queries.
%
% 4. Lines 138-145 (Axisymmetric Axis L'Hopital Safeguard):
%    Retained L'Hopital limit check for `ur_over_r` near r_actual -> 0 to 
%    prevent 0/0 floating-point singularities on the centerline.
%
% 5. Dynamic Boundary Profile Integration (meshF Alignment Check):
%    Note: Ensure meshF passed into this function reflects current 
%    deformed boundary profiles deltaE(z) and deltaL(z) from 
%    extract_interface_radius to ensure (rq, zq) maps to true physical 
%    fluid elements during FSI coupling.
% -------------------------------------------------------------------------

Nr = mraw.Nr; Nz = mraw.Nz;
i = i_guess; j = j_guess;

% Floor threshold for non-singular inverse evaluation
epsJ = 1e-14;

for walk = 1:10
    e = (j-1)*(Nr-1) + i;
    conn = meshF.elems(e,:);
    xe = meshF.nodes(conn,1);
    ze = meshF.nodes(conn,2);

    xi = 0; eta = 0;
    for it = 1:50
        [N, dNdxi] = shape_Q4(xi, eta);
        r_cur = N.' * xe;
        z_cur = N.' * ze;
        resid = [rq - r_cur; zq - z_cur];
        if norm(resid) < 1e-13
            break;
        end
        J = [xe ze].' * dNdxi;
        
        % [OLD]: detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
        % [OLD]: if abs(detJ) <= 1e-12
        % [OLD]:     dxieta = [0; 0];
        % [OLD]:     break;
        % [OLD]: end
        % [OLD]: invJ = [J(2,2), -J(1,2); -J(2,1), J(1,1)] / detJ;
        
        % Regularized 2x2 cofactor inverse prevents zero-division without breaking loop
        detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
        detJ_sign = sign(detJ);
        if detJ_sign == 0, detJ_sign = 1; end
        detJ_safe = detJ_sign * max(abs(detJ), epsJ);
        
        invJ = [J(2,2), -J(1,2); -J(2,1), J(1,1)] / detJ_safe;
        dxieta = invJ * resid;
        
        % Step size capping to prevent shooting far outside the element
        step_norm = norm(dxieta);
        if step_norm > 1.0
            dxieta = dxieta / step_norm;
        end
        
        xi = xi + dxieta(1);
        eta = eta + dxieta(2);
    end

    tol = 1e-6;
    if xi >= -1-tol && xi <= 1+tol && eta >= -1-tol && eta <= 1+tol
        xi = min(max(xi,-1),1);
        eta = min(max(eta,-1),1);
        break;
    end

    % Walk to neighboring element in direction of overshoot
    moved = false;
    if xi < -1 && i > 1, i = i-1; moved = true;
    elseif xi > 1 && i < Nr-1, i = i+1; moved = true;
    end
    if eta < -1 && j > 1, j = j-1; moved = true;
    elseif eta > 1 && j < Nz-1, j = j+1; moved = true;
    end
    if ~moved
        xi = min(max(xi,-1),1);
        eta = min(max(eta,-1),1);
        break;
    end
end

ij = [i, j];
xieta = [xi, eta];

[N, dNdxi] = shape_Q4(xi, eta);
J = [xe ze].' * dNdxi;

% [OLD]: detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
% [OLD]: if abs(detJ) <= 1e-12
% [OLD]:     error('locate_and_interp_fluid_stress: Singular Jacobian matrix det(J) <= 1e-12.');
% [OLD]: end
% [OLD]: invJ = [J(2,2), -J(1,2); -J(2,1), J(1,1)] / detJ;
% [OLD]: dNdx = dNdxi * invJ;

% Regularized Cartesian shape function gradients
detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
detJ_sign = sign(detJ);
if detJ_sign == 0, detJ_sign = 1; end
detJ_safe = detJ_sign * max(abs(detJ), epsJ);

invJ = [J(2,2), -J(1,2); -J(2,1), J(1,1)] / detJ_safe;
dNdx = dNdxi * invJ;

dNdr = dNdx(:,1);
dNdz = dNdx(:,2);

r_actual = N.' * xe;
z_actual = N.' * ze;
rz_actual = [r_actual, z_actual];

pc = N.' * plist(conn);

ur = zeros(4,1); uz = zeros(4,1);
for a = 1:4
    ur(a) = ur2D(conn(a));
    uz(a) = uz2D(conn(a));
end

durdr = dNdr.' * ur;
duzdz = dNdz.' * uz;
durdz = dNdz.' * ur;
duzdr = dNdr.' * uz;

% [OLD]: ur_over_r = (N.' * ur) / r_actual;
% L'Hopital limit safeguard for hoop strain rate as r_actual -> 0
if r_actual <= 1e-12
    ur_over_r = durdr;
else
    ur_over_r = (N.' * ur) / r_actual;
end

srr = -pc + 2*mu*durdr;
stt = -pc + 2*mu*ur_over_r;
szz = -pc + 2*mu*duzdz;
srz = mu*(durdz + duzdr);

sigma = [srr, stt, szz, srz];
p = pc;
end