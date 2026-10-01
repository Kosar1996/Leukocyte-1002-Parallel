function [Rsolid, Rfluid] = monolithic_analytical_residual_only( ...
    u, p, old, meshE, interfaceE, z, par, ...
    free, fixDofs, fixVals)

    ndof = size(meshE.nodes,1)*2;
    N = numel(z);

    u(fixDofs) = fixVals;
    p(1) = par.pIn;
    p(end) = par.pOut;

    [deltaE, UwE] = monolithic_interface_kinematics( ...
        meshE, u, old.uE, interfaceE, z, par);

    deltaL = old.deltaL;
    h = deltaE - deltaL;

    if any(h <= par.minGap)
        error('Gap violates minGap.');
    end

    UwL = par.UwL * ones(N,1);

    [Q, ~, ~, ~, tauE, ~, ~] = ...
        local_flux_and_shear(z, p, deltaL, deltaE, UwL, UwE, par);

    trE.normal = p;
    trE.tangent = -tauE;
    if use_exact_interface_in_monolithic(par)
        trE.z = z;
    end

    Fext = zeros(ndof,1);
    [Fext, ~] = apply_interface_traction(meshE, u, Fext, interfaceE, trE);

    Fint = assemble_finite_def_internal_force_only(meshE, u, par) + ...
        assemble_axisym_kelvin_voigt_viscous_force_only(meshE, u, old.uE, par);

    RsolidFull = Fint - Fext;
    Rsolid = RsolidFull(free);

    re    = deltaE(:);
    rl    = deltaL(:);
    reOld = old.deltaE(:);
    rlOld = old.deltaL(:);

    A = 0.5*(re.^2 - rl.^2);
    Aold = 0.5*(reOld.^2 - rlOld.^2);

    Ssrc = -par.SsrcFactor*(A - Aold) / par.dt;

    Rfluid = zeros(N-2,1);

    for i = 2:N-1
        Rfluid(i-1) = (Q(i) - Q(i-1))/par.dz - Ssrc(i);
    end
end
