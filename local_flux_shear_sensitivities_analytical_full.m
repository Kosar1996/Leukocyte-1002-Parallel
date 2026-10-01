function fs = local_flux_shear_sensitivities_analytical_full(z, p, rl, re, UwL, UwE, par)
% Analytical local sensitivities for the annular Navier-slip fluid map.
%
% Returns sparse global matrices for
%   Q(1:N-1), tauL(1:N), tauE(1:N)
% with respect to p, rl, re, UwL, and UwE.
%
% Convention:
%   rl = a = inner/leukocyte wall radius
%   re = b = outer/endothelium wall radius

    N  = numel(z);
    dzFace = diff(z(:));
    mu = par.mu;
    ll = par.slipL;
    le = par.slipE;

    dQdp    = sparse(N-1,N);
    dQdrl   = sparse(N-1,N);
    dQdre   = sparse(N-1,N);
    dQdUwL  = sparse(N-1,N);
    dQdUwE  = sparse(N-1,N);

    dtauLdp   = sparse(N,N);
    dtauLdrl  = sparse(N,N);
    dtauLdre  = sparse(N,N);
    dtauLdUwL = sparse(N,N);
    dtauLdUwE = sparse(N,N);

    dtauEdp   = sparse(N,N);
    dtauEdrl  = sparse(N,N);
    dtauEdre  = sparse(N,N);
    dtauEdUwL = sparse(N,N);
    dtauEdUwE = sparse(N,N);

    % Face flux derivatives. Each face depends on two neighboring wall radii,
    % two neighboring wall velocities, and the two adjacent pressures.
    for i = 1:N-1
        a  = 0.5*(rl(i)   + rl(i+1));
        b  = 0.5*(re(i)   + re(i+1));
        UL = 0.5*(UwL(i)  + UwL(i+1));
        UE = 0.5*(UwE(i)  + UwE(i+1));
        dz = dzFace(i);
        G  = (p(i+1)-p(i))/(mu*dz);

        d = annulus_local_derivatives(a, b, G, UL, UE, ll, le, mu);

        dQdp(i,i)   = dQdp(i,i)   - d.dQdG/(mu*dz);
        dQdp(i,i+1) = dQdp(i,i+1) + d.dQdG/(mu*dz);

        dQdrl(i,i)   = dQdrl(i,i)   + 0.5*d.dQda;
        dQdrl(i,i+1) = dQdrl(i,i+1) + 0.5*d.dQda;
        dQdre(i,i)   = dQdre(i,i)   + 0.5*d.dQdb;
        dQdre(i,i+1) = dQdre(i,i+1) + 0.5*d.dQdb;

        dQdUwL(i,i)   = dQdUwL(i,i)   + 0.5*d.dQdUL;
        dQdUwL(i,i+1) = dQdUwL(i,i+1) + 0.5*d.dQdUL;
        dQdUwE(i,i)   = dQdUwE(i,i)   + 0.5*d.dQdUE;
        dQdUwE(i,i+1) = dQdUwE(i,i+1) + 0.5*d.dQdUE;
    end

    % Nodal wall-shear derivatives. The pressure-gradient stencil here must
    % match local_flux_and_shear exactly.
    for i = 1:N
        a  = rl(i);
        b  = re(i);
        UL = UwL(i);
        UE = UwE(i);

        [G, pCols, pCoef] = global_1d_node_gradient_sensitivity(z, p, i, mu);

        d = annulus_local_derivatives(a, b, G, UL, UE, ll, le, mu);

        for k = 1:numel(pCols)
            j = pCols(k);
            dtauLdp(i,j) = dtauLdp(i,j) + d.dtauLdG*pCoef(k);
            dtauEdp(i,j) = dtauEdp(i,j) + d.dtauEdG*pCoef(k);
        end

        dtauLdrl(i,i)  = d.dtauLda;
        dtauLdre(i,i)  = d.dtauLdb;
        dtauLdUwL(i,i) = d.dtauLdUL;
        dtauLdUwE(i,i) = d.dtauLdUE;

        dtauEdrl(i,i)  = d.dtauEda;
        dtauEdre(i,i)  = d.dtauEdb;
        dtauEdUwL(i,i) = d.dtauEdUL;
        dtauEdUwE(i,i) = d.dtauEdUE;
    end

    fs.dQdp = dQdp;
    fs.dQdrl = dQdrl;
    fs.dQdre = dQdre;
    fs.dQdUwL = dQdUwL;
    fs.dQdUwE = dQdUwE;

    fs.dtauLdp = dtauLdp;
    fs.dtauLdrl = dtauLdrl;
    fs.dtauLdre = dtauLdre;
    fs.dtauLdUwL = dtauLdUwL;
    fs.dtauLdUwE = dtauLdUwE;

    fs.dtauEdp = dtauEdp;
    fs.dtauEdrl = dtauEdrl;
    fs.dtauEdre = dtauEdre;
    fs.dtauEdUwL = dtauEdUwL;
    fs.dtauEdUwE = dtauEdUwE;
end

function d = annulus_local_derivatives(a, b, G, UL, UE, ll, le, mu)
% Local analytical derivatives of the annular Stokes solution.
% u(r) = G*r^2/4 + C1*log(r) + C2, where G=(dp/dz)/mu.

    if a <= 0 || b <= 0 || b <= a
        error('Invalid annulus radii in analytical local sensitivity.');
    end

    m1 = log(a) + ll/a;
    m2 = log(b) + le/b;
    A1 = a^2/4 + ll*a/2;
    B1 = b^2/4 + le*b/2;

    M = [m1, 1; m2, 1];
    rhs = [UL - G*A1; UE - G*B1];
    C = M \ rhs;
    C1 = C(1);
    C2 = C(2); %#ok<NASGU>

    Ilog = 0.5*(b^2*log(b) - a^2*log(a)) - 0.25*(b^2 - a^2);
    I1   = 0.5*(b^2 - a^2);

    dIlog_da = -a*log(a);
    dIlog_db =  b*log(b);
    dI1_da   = -a;
    dI1_db   =  b;

    baseG = (b^4 - a^4)/16;

    % Derivatives of C=[C1;C2].
    dCdG  = M \ [-A1; -B1];
    dCdUL = M \ [1; 0];
    dCdUE = M \ [0; 1];

    dm1_da = 1/a - ll/a^2;
    dm2_db = 1/b - le/b^2;
    dA1_da = a/2 + ll/2;
    dB1_db = b/2 + le/2;

    dMda = [dm1_da, 0; 0, 0];
    dMdb = [0, 0; dm2_db, 0];
    drhs_da = [-G*dA1_da; 0];
    drhs_db = [0; -G*dB1_db];

    dCda = M \ (drhs_da - dMda*C);
    dCdb = M \ (drhs_db - dMdb*C);

    % Flux derivatives.
    dQdG  = baseG + dCdG(1)*Ilog + dCdG(2)*I1;
    dQdUL = dCdUL(1)*Ilog + dCdUL(2)*I1;
    dQdUE = dCdUE(1)*Ilog + dCdUE(2)*I1;
    dQda  = G*(-a^3/4) + dCda(1)*Ilog + C1*dIlog_da + dCda(2)*I1 + C(2)*dI1_da;
    dQdb  = G*( b^3/4) + dCdb(1)*Ilog + C1*dIlog_db + dCdb(2)*I1 + C(2)*dI1_db;

    % Shear derivatives.
    % tauL = mu*(G*a/2 + C1/a)
    % tauE = mu*(G*b/2 + C1/b)
    dtauLdG  = mu*(a/2 + dCdG(1)/a);
    dtauLdUL = mu*(dCdUL(1)/a);
    dtauLdUE = mu*(dCdUE(1)/a);
    dtauLda  = mu*(G/2 + dCda(1)/a - C1/a^2);
    dtauLdb  = mu*(dCdb(1)/a);

    dtauEdG  = mu*(b/2 + dCdG(1)/b);
    dtauEdUL = mu*(dCdUL(1)/b);
    dtauEdUE = mu*(dCdUE(1)/b);
    dtauEda  = mu*(dCda(1)/b);
    dtauEdb  = mu*(G/2 + dCdb(1)/b - C1/b^2);

    d = struct();
    d.dQdG = dQdG;
    d.dQda = dQda;
    d.dQdb = dQdb;
    d.dQdUL = dQdUL;
    d.dQdUE = dQdUE;

    d.dtauLdG = dtauLdG;
    d.dtauLda = dtauLda;
    d.dtauLdb = dtauLdb;
    d.dtauLdUL = dtauLdUL;
    d.dtauLdUE = dtauLdUE;

    d.dtauEdG = dtauEdG;
    d.dtauEda = dtauEda;
    d.dtauEdb = dtauEdb;
    d.dtauEdUL = dtauEdUL;
    d.dtauEdUE = dtauEdUE;
end

function [G, pCols, pCoef] = global_1d_node_gradient_sensitivity(z, p, i, mu)
    z = z(:);
    p = p(:);
    N = numel(z);
    if i <= 1
        denom = z(2)-z(1);
        pCols = [1, 2];
        pCoef = [-1/(mu*denom), 1/(mu*denom)];
    elseif i >= N
        denom = z(N)-z(N-1);
        pCols = [N-1, N];
        pCoef = [-1/(mu*denom), 1/(mu*denom)];
    else
        denom = z(i+1)-z(i-1);
        pCols = [i-1, i+1];
        pCoef = [-1/(mu*denom), 1/(mu*denom)];
    end
    G = p(pCols(1))*pCoef(1) + p(pCols(2))*pCoef(2);
end