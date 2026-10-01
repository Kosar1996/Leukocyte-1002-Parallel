function dF = local_dF_from_dof(alpha, N, dNdX, Rg)
% Returns the exact variation dF associated with one local DOF alpha.
%
% Local DOF numbering:
%   1 = ur1, 2 = uz1, 3 = ur2, 4 = uz2, 5 = ur3, 6 = uz3, 7 = ur4, 8 = uz4

    aNode = ceil(alpha/2);
    isRadial = mod(alpha,2)==1;

    dF = zeros(3,3);

    if isRadial
        % Only r-components vary
        d_r_g   = N(aNode);
        d_drdR  = dNdX(aNode,1);
        d_drdZ  = dNdX(aNode,2);

        % dF(1,1) = d_drdR;
        % dF(1,3) = d_drdZ;
        % dF(2,2) = d_r_g / Rg;
        % ========================= NEW CODE ==============================
        epsR = 1e-14;
        Rg_eff = max(Rg, epsR);

        dF(1,1) = d_drdR;
        dF(1,3) = d_drdZ;

        % L'Hopital Limit for dF(2,2) variation (lim_{R->0} Na/R = dNa/dR)
        if Rg < 1e-10
            dF(2,2) = d_drdR;
        else
            dF(2,2) = d_r_g / Rg_eff;
        end
        % =================================================================

    else
        % Only z-components vary
        d_dzdR = dNdX(aNode,1);
        d_dzdZ = dNdX(aNode,2);

        dF(3,1) = d_dzdR;
        dF(3,3) = d_dzdZ;
    end
end