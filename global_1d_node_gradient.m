function dpdz = global_1d_node_gradient(z, p, i)
    z = z(:);
    p = p(:);
    N = numel(z);
    if i <= 1
        dpdz = (p(2)-p(1))/(z(2)-z(1));
    elseif i >= N
        dpdz = (p(N)-p(N-1))/(z(N)-z(N-1));
    else
        dpdz = (p(i+1)-p(i-1))/(z(i+1)-z(i-1));
    end
end