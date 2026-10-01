function plot_select_native2d_stress(out,plotstep)
if ~isfield(out.fluidHist{plotstep}, 'P') || isempty(out.fluidHist{plotstep}.P) || ...
        ~any(isfinite(out.fluidHist{plotstep}.P(:)))
    warning('No finite native 2D pressure field is available to plot.');
    return;
end

figure;
ax = gca;
set(ax, 'FontSize', 22);
%draw fluid stress/pressure
if is_hybrid_fluid_plot(out)
    contour_masked_regular_field_select(ax, out, plotstep, 1.0, 48);
end
colorbar;
hold(ax, 'on');
% draw solid stress/pressure
stressE = recover_nodal_stress_axisym(out.meshE, out.stateHist{plotstep}.uE, out.par);
plot_select_nodal_stress_contour(out.meshE, out.stateHist{plotstep}.uE, -(stressE.sigma_rr+stressE.sigma_zz+stressE.sigma_tt)/3);
stressL = recover_nodal_stress_axisym(out.meshL, out.stateHist{plotstep}.uL, out.par);
plot_select_nodal_stress_contour(out.meshL, out.stateHist{plotstep}.uL, -(stressL.sigma_rr+stressL.sigma_zz+stressL.sigma_tt)/3);
draw_select_fluid_boundaries(ax, out,plotstep);

xlabel(ax, 'r [\mum]');
ylabel(ax, 'z [\mum]');
title(ax, sprintf('Hybrid pressure P(r,z), t = %.4g s', ...
    out.t(plotstep)));
axis(ax, 'equal');
axis(ax, 'tight');
xlim([0 4]);
box(ax, 'on');
end

function contour_masked_regular_field_select(ax, out, plotstep, valueScale, nLevels)
% Plot the hybrid field on a regular r-z image grid and mask the deformed
% solids. This avoids drawing artificial curvilinear cells that connect the
% 1D gap strip directly to the 2D exterior reservoir across solid caps.

%this version uses the pressures interpolated from nodal velocity field
native2D.P = reshape(out.fluidHist{plotstep}.pCellNode,out.fluidHist{plotstep}.meshF.Nr-1,out.fluidHist{plotstep}.meshF.Nz-1); 
native2D.R = reshape(out.fluidHist{plotstep}.centerNode(:,1),out.fluidHist{plotstep}.meshF.Nr-1,out.fluidHist{plotstep}.meshF.Nz-1);
native2D.Z = reshape(out.fluidHist{plotstep}.centerNode(:,2),out.fluidHist{plotstep}.meshF.Nr-1,out.fluidHist{plotstep}.meshF.Nz-1);

%this version uses pressure already stored in each node in the original
%code and would recover identical result as the original code
% native2D.P = reshape(out.fluidHist{plotstep}.pCell,out.fluidHist{plotstep}.meshF.Nr,out.fluidHist{plotstep}.meshF.Nz);
% native2D.R = reshape(out.fluidHist{plotstep}.cellCenter(:,1),out.fluidHist{plotstep}.meshF.Nr,out.fluidHist{plotstep}.meshF.Nz);
% native2D.Z = reshape(out.fluidHist{plotstep}.cellCenter(:,2),out.fluidHist{plotstep}.meshF.Nr,out.fluidHist{plotstep}.meshF.Nz);

valid = isfinite(native2D.R) & isfinite(native2D.Z) & isfinite(native2D.P);
if nnz(valid) < 3
    warning('Not enough finite hybrid field points to plot.');
    return;
end

rVals = native2D.R(valid);
zVals = native2D.Z(valid);
fVals = valueScale * native2D.P(valid);

nrPlot = 260;
nzPlot = 260;
rGrid = linspace(min(rVals), max(rVals), nrPlot);
zGrid = linspace(min(zVals), max(zVals), nzPlot);
[RGrid, ZGrid] = meshgrid(rGrid, zGrid);

F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
FGrid = F(RGrid, ZGrid);

statePlot = state_for_plot_at_step(out, plotstep);
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


