function [Q, dQdpL, dQdpR, tauL, tauE, uzL, uzE] = local_flux_and_shear(z, p, rl, re, UwL, UwE, par)
    % Navier-slip solution for velocity, flux, and stress
    
    N = par.NzFluid;
    dzFace = diff(z(:));
    ll = par.slipL;
    le = par.slipE;
    mu = par.mu;
    
    Q = zeros(N-1,1); % axial flux through i
    dQdpL = zeros(N-1,1); % derivative of that face flux with respect to left pressure
    dQdpR = zeros(N-1,1); % derivative with respect to right pressure
    tauL = zeros(N,1); % axial shear stress at leukocyte wall
    tauE = zeros(N,1); % axial shear stress at endothelium wall
    uzL = zeros(N,1);  % axial fluid velocity at the leukocyte surface
    uzE = zeros(N,1);  % axial fluid velocity at the endothelium surface

    if isfield(par, 'noLeukocyte') && par.noLeukocyte
        for i = 1:N-1
            b = 0.5*(re(i)+re(i+1));
            dz = dzFace(i);
            dpdz = (p(i+1)-p(i))/dz;
            G = dpdz/mu;
            UwEi = 0.5*(UwE(i) + UwE(i+1));

            dQdG = -b^4/16 - le*b^3/4;
            Q(i) = dQdG*G + UwEi*b^2/2;
            dQdpL(i) = -dQdG/(mu*dz);
            dQdpR(i) =  dQdG/(mu*dz);
        end

        for i = 1:N
            b = re(i);
            dpdz = global_1d_node_gradient(z, p, i);
            G = dpdz/mu;
            uzL(i) = UwE(i) - G*b^2/4 - le*G*b/2;
            uzE(i) = UwE(i);
            tauL(i) = 0;
            tauE(i) = mu*G*b/2;
        end
        return;
    end
    
    for i = 1:N-1
        a = 0.5*(rl(i)+rl(i+1));
        b = 0.5*(re(i)+re(i+1));
        dz = dzFace(i);
        dpdz = (p(i+1)-p(i))/dz; 
        G = dpdz/mu;
        % Solve u_z(r) = (dp/dz/4) r^2 + C1 log r + C2
        % with Navier-slip BCs:
        % u(a) - UwL = - l_l u_r(a)
        % u(b)     = - l_e u_r(b)

        UwLi = 0.5*(UwL(i) + UwL(i+1));
        UwEi = 0.5*(UwE(i) + UwE(i+1));
    
        rhs = [UwLi - (G/4)*a^2 - ll*(G*a/2);
               UwEi        - (G/4)*b^2 - le*(G*b/2)];
        M = [log(a)+ll/a, 1; log(b)+le/b, 1];
        C = M\rhs;

        C1 = C(1); 
        C2 = C(2);
    
        % Flux Q = int_a^b r*u(r) dr
        Q(i) = integral_flux_annulus(a,b,G,C1,C2);
    
        drhsdG = [-(a^2/4 + ll*a/2); 
                  -(b^2/4 + le*b/2)];
        dCdG = M \ drhsdG;
        dC1dG=dCdG(1);
        dC2dG=dCdG(2);
        
        Ilog = 0.5*(b^2*log(b) - a^2*log(a)) - 0.25*(b^2 - a^2);
        I1   = 0.5*(b^2 - a^2);
        
        dQdG = (b^4 - a^4)/16 + dC1dG*Ilog + dC2dG*I1;
        dQdpL(i) = -dQdG/(mu*dz);
        dQdpR(i) =  dQdG/(mu*dz);
    end
    
    % Shear and interfacial axial speed
    for i = 1:N
        a = rl(i);
        b = re(i);
        dpdz = global_1d_node_gradient(z, p, i);
        UwLi = UwL(i);
        UwEi = UwE(i);
        G = dpdz/mu;
        M = [log(a)+ll/a, 1; log(b)+le/b, 1];
        rhs = [UwLi - (G/4)*a^2 - ll*(G*a/2);
               UwEi - (G/4)*b^2 - le*(G*b/2)];
        C = M\rhs;
        C1 = C(1); C2 = C(2);
    uzL(i) = (G/4)*a^2 + C1*log(a) + C2;
    uzE(i) = (G/4)*b^2 + C1*log(b) + C2;
    tauL(i) = mu*(G*a/2 + C1/a);
    tauE(i) = mu*(G*b/2 + C1/b);
    end
end