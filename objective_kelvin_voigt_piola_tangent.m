function dPvisc = objective_kelvin_voigt_piola_tangent(dF, data, par)
    I3 = eye(3);
    dFdot = dF / par.dt;
    dL = dFdot * data.Finv - data.L * dF * data.Finv;
    dD = 0.5 * (dL + dL.');
    dtrD = dD(1,1) + dD(2,2) + dD(3,3);
    ddevD = dD - (dtrD/3) * I3;
    dTvisc = 2 * par.etaE * ddevD + data.etaBulk * dtrD * I3;

    trFinvDF = sum(sum(data.Finv.' .* dF));
    dJ = data.J * trFinvDF;
    dFinvT = -data.FinvT * dF.' * data.FinvT;

    dPvisc = dJ * data.Tvisc * data.FinvT + ...
        data.J * dTvisc * data.FinvT + ...
        data.J * data.Tvisc * dFinvT;
end
