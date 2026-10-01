function [solidTarget, fluidTarget] = monolithic_residual_targets(par)
    solidTarget = par.solidAbsTol;
    if isfield(par, 'monoSolidAbsTol')
        solidTarget = par.monoSolidAbsTol;
    end

    fluidTarget = par.tolNewtonFluid;
    if isfield(par, 'monoFluidAbsTol')
        fluidTarget = par.monoFluidAbsTol;
    end
end