function sens = local_flux_shear_sensitivities(z, p, rl, re, UwL, UwE, par)

    N = par.NzFluid;
    dz = par.dz;
    mu = par.mu;
    ll = par.slipL;
    le = par.slipE;

    Q = zeros(N-1,1);

    dQdpL  = zeros(N-1,1);
    dQdpR  = zeros(N-1,1);
    dQdb   = zeros(N-1,1);
    dQdUwE = zeros(N-1,1);

    tauE = zeros(N,1);

    dtauEdp    = sparse(N,N);
    dtauEdb    = sparse(N,N);
    dtauEdUwE  = sparse(N,N);

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        for i = 1:N-1
            b = 0.5*(re(i)+re(i+1));
            dpdz = (p(i+1)-p(i))/dz;
            G = dpdz/mu;
            UwEi = 0.5*(UwE(i)+UwE(i+1));

            dQdG = -b^4/16 - le*b^3/4;
            Q(i) = dQdG*G + UwEi*b^2/2;

            dQdpL(i) = -dQdG/(mu*dz);
            dQdpR(i) =  dQdG/(mu*dz);
            dQdb(i) = (-G*b^3/4 - 3*le*G*b^2/4) + UwEi*b;
            dQdUwE(i) = b^2/2;
        end

        for i = 1:N
            if i == 1
                b = re(i);
                dpdz = (p(i+1)-p(i))/dz;
                p_ids = [i, i+1];
                dp_coeff = [-1/dz, 1/dz];
            elseif i == N
                b = re(i);
                dpdz = (p(i)-p(i-1))/dz;
                p_ids = [i-1, i];
                dp_coeff = [-1/dz, 1/dz];
            else
                b = re(i);
                dpdz = (p(i+1)-p(i-1))/(2*dz);
                p_ids = [i-1, i+1];
                dp_coeff = [-1/(2*dz), 1/(2*dz)];
            end

            G = dpdz/mu;
            tauE(i) = mu*G*b/2;

            for k = 1:numel(p_ids)
                dtauEdp(i,p_ids(k)) = dtauEdp(i,p_ids(k)) + (b/2) * dp_coeff(k);
            end
            dtauEdb(i,i) = mu*G/2;
        end

        sens.Q = Q;
        sens.dQdpL = dQdpL;
        sens.dQdpR = dQdpR;
        sens.dQdb = dQdb;
        sens.dQdUwE = dQdUwE;

        sens.tauE = tauE;
        sens.dtauEdp = dtauEdp;
        sens.dtauEdb = dtauEdb;
        sens.dtauEdUwE = dtauEdUwE;
        return;
    end

    for i = 1:N-1

        a = 0.5*(rl(i)+rl(i+1));
        b = 0.5*(re(i)+re(i+1));

        dpdz = (p(i+1)-p(i))/dz;
        G = dpdz/mu;

        UwLi = 0.5*(UwL(i)+UwL(i+1));
        UwEi = 0.5*(UwE(i)+UwE(i+1));

        M = [log(a)+ll/a, 1;
             log(b)+le/b, 1];

        rhs = [UwLi - (G/4)*a^2 - ll*(G*a/2);
               UwEi - (G/4)*b^2 - le*(G*b/2)];

        C = M\rhs;

        C1 = C(1);
        C2 = C(2);

        Q(i) = integral_flux_annulus(a,b,G,C1,C2);

        Ilog = 0.5*(b^2*log(b)-a^2*log(a)) - 0.25*(b^2-a^2);
        I1   = 0.5*(b^2-a^2);

        % dQ/dG
        drhsdG = [-(a^2/4 + ll*a/2);
                  -(b^2/4 + le*b/2)];

        dCdG = M \ drhsdG;

        dQdG = (b^4-a^4)/16 + dCdG(1)*Ilog + dCdG(2)*I1;

        dQdpL(i) = -dQdG/(mu*dz);
        dQdpR(i) =  dQdG/(mu*dz);

        % dQ/dUwE
        dCdUwE = M \ [0;1];

        dQdUwE(i) = dCdUwE(1)*Ilog + dCdUwE(2)*I1;

        % dQ/db
        dMdb = [0, 0;
                1/b - le/b^2, 0];

        drhsdb = [0;
                 -(G/2)*b - le*(G/2)];

        dCdb = M \ (drhsdb - dMdb*C);

        dIlogdb = b*log(b);
        dI1db = b;

        dQdb(i) = ...
            (G/4)*b^3 ...
            + dCdb(1)*Ilog + C1*dIlogdb ...
            + dCdb(2)*I1   + C2*dI1db;
    end

    for i = 1:N

        if i == 1
            a = rl(i);
            b = re(i);
            dpdz = (p(i+1)-p(i))/dz;
            p_ids = [i, i+1];
            dp_coeff = [-1/dz, 1/dz];

        elseif i == N
            a = rl(i);
            b = re(i);
            dpdz = (p(i)-p(i-1))/dz;
            p_ids = [i-1, i];
            dp_coeff = [-1/dz, 1/dz];

        else
            a = rl(i);
            b = re(i);
            dpdz = (p(i+1)-p(i-1))/(2*dz);
            p_ids = [i-1, i+1];
            dp_coeff = [-1/(2*dz), 1/(2*dz)];
        end

        G = dpdz/mu;

        UwLi = UwL(i);
        UwEi = UwE(i);

        M = [log(a)+ll/a, 1;
             log(b)+le/b, 1];

        rhs = [UwLi - (G/4)*a^2 - ll*(G*a/2);
               UwEi - (G/4)*b^2 - le*(G*b/2)];

        C = M\rhs;
        C1 = C(1);

        tauE(i) = mu*(G*b/2 + C1/b);

        % dtauE/dG
        drhsdG = [-(a^2/4 + ll*a/2);
                  -(b^2/4 + le*b/2)];

        dCdG = M \ drhsdG;

        dtauEdG = mu*(b/2 + dCdG(1)/b);

        for k = 1:numel(p_ids)
            dtauEdp(i,p_ids(k)) = dtauEdp(i,p_ids(k)) + dtauEdG * dp_coeff(k) / mu;
        end

        % dtauE/dUwE
        dCdUwE = M \ [0;1];

        dtauEdUwE(i,i) = mu*dCdUwE(1)/b;

        % dtauE/db
        dMdb = [0, 0;
                1/b - le/b^2, 0];

        drhsdb = [0;
                 -(G/2)*b - le*(G/2)];

        dCdb = M \ (drhsdb - dMdb*C);

        dtauEdb(i,i) = mu*(G/2 + dCdb(1)/b - C1/b^2);
    end

    sens.Q = Q;
    sens.dQdpL = dQdpL;
    sens.dQdpR = dQdpR;
    sens.dQdb = dQdb;
    sens.dQdUwE = dQdUwE;

    sens.tauE = tauE;
    sens.dtauEdp = dtauEdp;
    sens.dtauEdb = dtauEdb;
    sens.dtauEdUwE = dtauEdUwE;
end