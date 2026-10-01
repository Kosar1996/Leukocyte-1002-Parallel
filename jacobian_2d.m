function [Jm, dNdx, detJ] = jacobian_2d(Xe, dNdxi)
    Jm = Xe.' * dNdxi;
    detJ = det(Jm);
    if detJ <= 0
        error('Non-positive element Jacobian.');
    end
    dNdx = dNdxi / Jm;
end