function tf = use_global2d_pressure_traction(par)
    tf = isfield(par, 'useGlobal2DPressureTraction') && ...
        par.useGlobal2DPressureTraction;
end
