function fs = local_flux_shear_sensitivities_fd_full(z, p, rl, re, UwL, UwE, par, Q0, tauL0, tauE0)
% Cheap local/global finite-difference sensitivities of the fluid-only map.
% This finite-differences only the 1D Reynolds/annular solution. It does NOT
% reassemble either solid, so it is much cheaper than letting fsolve finite-
% difference the full two-solid residual.

    N = numel(z);
    if nargin < 8 || isempty(Q0)
        [Q0, ~, ~, tauL0, tauE0, ~, ~] = local_flux_and_shear(z, p, rl, re, UwL, UwE, par);
    end

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

    epsP = 1e-6 * max(1, max(abs(p)));
    epsR = 1e-8 * max([max(abs(rl)), max(abs(re)), 1e-6]);
    epsU = 1e-8 * max([max(abs(UwL)), max(abs(UwE)), 1e-6]);

    % Pressure perturbations.
    for j = 1:N
        pp = p;
        pp(j) = pp(j) + epsP;
        [Qp, ~, ~, tauLp, tauEp, ~, ~] = local_flux_and_shear(z, pp, rl, re, UwL, UwE, par);
        dQdp(:,j) = (Qp - Q0) / epsP;
        dtauLdp(:,j) = (tauLp - tauL0) / epsP;
        dtauEdp(:,j) = (tauEp - tauE0) / epsP;
    end

    % Inner radius / leukocyte wall perturbations.
    for j = 1:N
        rlp = rl;
        rlp(j) = rlp(j) + epsR;
        if any(re - rlp <= par.minGap)
            rlp(j) = rl(j) - epsR;
            sgn = -1;
        else
            sgn = 1;
        end
        [Qp, ~, ~, tauLp, tauEp, ~, ~] = local_flux_and_shear(z, p, rlp, re, UwL, UwE, par);
        dQdrl(:,j) = sgn*(Qp - Q0) / epsR;
        dtauLdrl(:,j) = sgn*(tauLp - tauL0) / epsR;
        dtauEdrl(:,j) = sgn*(tauEp - tauE0) / epsR;
    end

    % Outer radius / endothelium wall perturbations.
    for j = 1:N
        rep = re;
        rep(j) = rep(j) + epsR;
        [Qp, ~, ~, tauLp, tauEp, ~, ~] = local_flux_and_shear(z, p, rl, rep, UwL, UwE, par);
        dQdre(:,j) = (Qp - Q0) / epsR;
        dtauLdre(:,j) = (tauLp - tauL0) / epsR;
        dtauEdre(:,j) = (tauEp - tauE0) / epsR;
    end

    % Inner wall axial velocity perturbations.
    for j = 1:N
        up = UwL;
        up(j) = up(j) + epsU;
        [Qp, ~, ~, tauLp, tauEp, ~, ~] = local_flux_and_shear(z, p, rl, re, up, UwE, par);
        dQdUwL(:,j) = (Qp - Q0) / epsU;
        dtauLdUwL(:,j) = (tauLp - tauL0) / epsU;
        dtauEdUwL(:,j) = (tauEp - tauE0) / epsU;
    end

    % Outer wall axial velocity perturbations.
    for j = 1:N
        up = UwE;
        up(j) = up(j) + epsU;
        [Qp, ~, ~, tauLp, tauEp, ~, ~] = local_flux_and_shear(z, p, rl, re, UwL, up, par);
        dQdUwE(:,j) = (Qp - Q0) / epsU;
        dtauLdUwE(:,j) = (tauLp - tauL0) / epsU;
        dtauEdUwE(:,j) = (tauEp - tauE0) / epsU;
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