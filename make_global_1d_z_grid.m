% Axial fluid grid (shared with interface interpolation locations)
function z = make_global_1d_z_grid(par)
    if ~(isfield(par, 'useGlobal1DPressure') && par.useGlobal1DPressure) || ...
            ~(isfield(par, 'useGlobal1DCoarseEdgeMesh') && par.useGlobal1DCoarseEdgeMesh)
        z = linspace(par.zMin, par.zMax, par.NzFluid).';
        return;
    end

    fineWindow = [-2e-6, 6e-6];
    if isfield(par, 'global1DFineWindow') && numel(par.global1DFineWindow) == 2
        fineWindow = sort(par.global1DFineWindow(:)).';
    end
    fineWindow(1) = max(fineWindow(1), par.zMin);
    fineWindow(2) = min(fineWindow(2), par.zMax);

    coarseDz = 1e-6;
    fineDz = 0.2e-6;
    if isfield(par, 'global1DCoarseDz') && isfinite(par.global1DCoarseDz) && par.global1DCoarseDz > 0
        coarseDz = par.global1DCoarseDz;
    end
    if isfield(par, 'global1DFineDz') && isfinite(par.global1DFineDz) && par.global1DFineDz > 0
        fineDz = par.global1DFineDz;
    end

    z = sort([ ...
        segment_grid(par.zMin, fineWindow(1), coarseDz)
        segment_grid(fineWindow(1), fineWindow(2), fineDz)
        segment_grid(fineWindow(2), par.zMax, coarseDz)
        par.zMin
        par.zMax
        fineWindow(:)
        0
        4e-6]);

    z = z(z >= par.zMin - 100*eps(max(abs(par.zMin),1)) & ...
          z <= par.zMax + 100*eps(max(abs(par.zMax),1)));
    z(abs(z) < 1e-15) = 0;
    z(abs(z - 4e-6) < 1e-15) = 4e-6;
    keep = [true; diff(z) > 1e-15];
    z = z(keep);
end

function z = segment_grid(a, b, dz)
    if b < a
        z = zeros(0,1);
        return;
    end
    if abs(b-a) < 1e-18
        z = a;
        return;
    end
    n = max(1, ceil((b-a)/dz));
    z = linspace(a, b, n+1).';
end