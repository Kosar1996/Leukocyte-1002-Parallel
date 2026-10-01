function deltaE = exact_endothelium_delta_for_gap(s, par)
    deltaE = [];
    if use_exact_endothelium_fluid_domain(par) && ...
            isfield(s, 'deltaEFluidBoundary')
        deltaE = s.deltaEFluidBoundary(:);
    elseif isfield(s, 'deltaESolid')
        deltaE = s.deltaESolid(:);
    end
end