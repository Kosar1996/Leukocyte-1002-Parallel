function [N, dNdxi, w] = q4_shape(xi, eta, gw)
    N = 0.25*[(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
    dNdxi = 0.25*[ -(1-eta), -(1-xi);
                    +(1-eta), -(1+xi);
                    +(1+eta), +(1+xi);
                    -(1+eta), +(1-xi) ];
    w = gw;
end
