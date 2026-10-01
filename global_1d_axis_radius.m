function r0 = global_1d_axis_radius(par)
    r0 = 1e-9;
    if isfield(par, 'global1DAxisRadius') && ...
            isfinite(par.global1DAxisRadius) && par.global1DAxisRadius > 0
        r0 = par.global1DAxisRadius;
    end
end