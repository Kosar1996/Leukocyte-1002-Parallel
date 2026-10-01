function tf = use_hybrid_gap1d_exterior2d_fluid(par)
    tf = isfield(par, 'useHybridGap1DExterior2DFluid') && ...
        par.useHybridGap1DExterior2DFluid;
end