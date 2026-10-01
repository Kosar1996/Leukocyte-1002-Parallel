function contour_masked_regular_field(ax, out, fieldNative, valueScale, nLevels)
% Plot the hybrid field on a regular r-z image grid and mask the deformed
% solids. This avoids drawing artificial curvilinear cells that connect the
% 1D gap strip directly to the 2D exterior reservoir across solid caps.
    native2D = out.native2D;
    valid = isfinite(native2D.R) & isfinite(native2D.Z) & isfinite(fieldNative);
    if nnz(valid) < 3
        warning('Not enough finite hybrid field points to plot.');
        return;
    end

    rVals = native2D.R(valid);
    zVals = native2D.Z(valid);
    fVals = valueScale * fieldNative(valid);

    nrPlot = 260;
    nzPlot = 260;
    rGrid = linspace(min(rVals), max(rVals), nrPlot);
    zGrid = linspace(min(zVals), max(zVals), nzPlot);
    [RGrid, ZGrid] = meshgrid(rGrid, zGrid);

    F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
    FGrid = F(RGrid, ZGrid);

    statePlot = final_state_for_plot(out);
    solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid);
    FGrid(solidMask) = NaN;

    finiteVals = FGrid(isfinite(FGrid));
    if isempty(finiteVals)
        warning('No finite hybrid field values remain after solid masking.');
        return;
    end
    if max(finiteVals) > min(finiteVals)
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, nLevels, 'LineColor', 'none');
    else
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 1, 'LineColor', 'none');
    end
end

function solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid)
    solidMask = false(size(RGrid));
    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL)
        solidMask = solidMask | local_deformed_solid_mask( ...
            out.meshL, statePlot.uL, RGrid, ZGrid);
    end
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        solidMask = solidMask | local_deformed_solid_mask( ...
            out.meshE, statePlot.uE, RGrid, ZGrid);
    end
end

function solidMask = local_deformed_solid_mask(mesh, u, RGrid, ZGrid)
    solidMask = false(size(RGrid));
    if isempty(mesh) || isempty(u) || ~isfield(mesh, 'nodes') || ...
            ~isfield(mesh, 'conn') || numel(u) < 2*size(mesh.nodes,1)
        return;
    end

    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    for e = 1:size(mesh.conn,1)
        ids = mesh.conn(e,:);
        rv = rDef(ids);
        zv = zDef(ids);
        inBox = RGrid >= min(rv) & RGrid <= max(rv) & ...
                ZGrid >= min(zv) & ZGrid <= max(zv);
        if any(inBox(:))
            localMask = false(size(solidMask));
            localMask(inBox) = inpolygon(RGrid(inBox), ZGrid(inBox), rv, zv);
            solidMask = solidMask | localMask;
        end
    end
end

