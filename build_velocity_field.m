function [R, Z, Uz] = build_velocity_field(z, p, rl, re, UwL, UwE, par)

    Nr = par.NrPlot;
    Nz = numel(z);

    R  = zeros(Nr, Nz);
    Z  = zeros(Nr, Nz);
    Uz = zeros(Nr, Nz);

    ll = par.slipL;
    le = par.slipE;
    mu = par.mu;
    dz = z(2)-z(1);

    for j = 1:Nz
        a = rl(j);
        b = re(j);

        if j == 1
            dpdz = (p(j+1)-p(j))/dz;
        elseif j == Nz
            dpdz = (p(j)-p(j-1))/dz;
        else
            dpdz = (p(j+1)-p(j-1))/(2*dz);
        end

        G = dpdz/mu;

        UwLj = UwL(j);
        UwEj = UwE(j);

        if isfield(par, 'noLeukocyte') && par.noLeukocyte
            r = linspace(0, b, Nr).';
            uz = (G/4)*(r.^2 - b^2) - le*G*b/2 + UwEj;
            R(:,j)  = r;
            Z(:,j)  = z(j);
            Uz(:,j) = uz;
            continue;
        end

        rhs = [UwLj - (G/4)*a^2 - ll*(G*a/2);
               UwEj - (G/4)*b^2 - le*(G*b/2)];
        M = [log(a)+ll/a, 1;
             log(b)+le/b, 1];
        C = M\rhs;
        C1 = C(1);
        C2 = C(2);

        r = linspace(a, b, Nr).';
        uz = (G/4)*r.^2 + C1*log(r) + C2;
        R(:,j)  = r;
        Z(:,j)  = z(j);
        Uz(:,j) = uz;
    end
end
