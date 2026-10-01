function [N, dNdxi] = shape_Q4(xi, eta)
    N = 0.25 * [
        (1-xi)*(1-eta);
        (1+xi)*(1-eta);
        (1+xi)*(1+eta);
        (1-xi)*(1+eta)];

    dN_dxi = 0.25 * [
        -(1-eta);
         (1-eta);
         (1+eta);
        -(1+eta)];

    dN_deta = 0.25 * [
        -(1-xi);
        -(1+xi);
         (1+xi);
         (1-xi)];

    dNdxi = [dN_dxi, dN_deta];
end