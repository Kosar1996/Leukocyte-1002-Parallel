function R = fluid_residual_only(Q, p, S, dz, par)
% how much current pressure field violates mass conservation + boundary conditions
    N = par.NzFluid;
    R = zeros(N,1);
    if isfield(par, 'zGrid') && numel(par.zGrid) == N
        dzFace = diff(par.zGrid(:));
        dzControl = global_1d_control_lengths(par.zGrid(:));
    else
        dzFace = dz * ones(N-1,1); %#ok<NASGU>
        dzControl = dz * ones(N,1);
    end
    % Pressure boundary conditions
    R(1) = p(1)-par.pIn; 
    R(end) = p(end)-par.pOut;
    % Mass conservation
    for i=2:N-1
        R(i) = (Q(i)-Q(i-1))/dzControl(i) - S(i);
    end
end