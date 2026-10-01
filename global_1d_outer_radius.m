function rMax = global_1d_outer_radius(par)
    rMax = 15e-6;
    if isfield(par, 'global1DOuterRadius') && ...
            isfinite(par.global1DOuterRadius) && par.global1DOuterRadius > 0
        rMax = par.global1DOuterRadius;
    end
end