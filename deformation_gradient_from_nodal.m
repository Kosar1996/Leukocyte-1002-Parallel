function F = deformation_gradient_from_nodal(~, rnod, znod, N, dNdX, Rg)
    rg = N * rnod;
    if rg <= 0
        error('Non-positive radius encountered in viscoelastic element.');
    end

    drdR = dNdX(:,1).' * rnod;
    drdZ = dNdX(:,2).' * rnod;
    dzdR = dNdX(:,1).' * znod;
    dzdZ = dNdX(:,2).' * znod;

    F = [drdR,   0,    drdZ;
           0,   rg/Rg, 0;
         dzdR,   0,    dzdZ];
end