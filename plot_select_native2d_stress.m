%% =========================================================================
% AUDITED & ROBUST NATIVE 2D STRESS PLOTTING FUNCTION
% Renders 2x2 Cauchy stress component tiles across solid and fluid domains.
% =========================================================================

function plot_select_native2d_stress(out, plotstep, varargin)
    if ~isfield(out, 'z') && isfield(out, 'meshE') && isfield(out.meshE, 'nodes')
        out.z = out.meshE.nodes(:,2);
    end

    if nargin < 2 || isempty(plotstep)
        plotstep = out.stopStep;
    end

    if plotstep < 1 || plotstep > numel(out.fluidHist)
        warning('plot_select_native2d_stress: plotstep %d is out of bounds.', plotstep);
        return;
    end

    fluid = out.fluidHist{plotstep};

    if isfield(fluid, 'sigmaCellNode') && ~isempty(fluid.sigmaCellNode)
        fluidSigma  = fluid.sigmaCellNode;
        fluidCenter = fluid.centerNode;
    elseif isfield(fluid, 'sigmaCell') && ~isempty(fluid.sigmaCell)
        fluidSigma  = fluid.sigmaCell;
        fluidCenter = fluid.cellCenter;
    else
        warning('plot_select_native2d_stress: No recovered fluid stress field found.');
        return;
    end

    if size(fluidSigma, 2) < 4 || ~any(isfinite(fluidSigma(:)))
        warning('plot_select_native2d_stress: Invalid fluid stress matrix.');
        return;
    end

    st = out.stateHist{plotstep};

    parL = out.par;
    if isfield(out.par, 'GL') && isfinite(out.par.GL), parL.Ge = out.par.GL; end
    if isfield(out.par, 'KL') && isfinite(out.par.KL), parL.Ke = out.par.KL; end
    if isfield(out.par, 'etaL') && isfinite(out.par.etaL), parL.etaE = out.par.etaL; end

    parE = out.par;
    if ~isfield(parE, 'etaE') && isfield(out.par, 'etaE'), parE.etaE = out.par.etaE; end

    hasPrevE = isfield(st, 'uEPrev') && ~isempty(st.uEPrev) && numel(st.uEPrev) == numel(st.uE);
    hasPrevL = isfield(st, 'uLPrev') && ~isempty(st.uLPrev) && numel(st.uLPrev) == numel(st.uL);
    hasDt = isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && ...
        isfinite(out.dtHist(plotstep)) && out.dtHist(plotstep) > 0;

    useViscoRecovery = hasPrevE && hasPrevL && hasDt;

    if useViscoRecovery
        dtStep = out.dtHist(plotstep);
        stressE = recover_nodal_stress_axisym_viscoelastic(out.meshE, st.uE, st.uEPrev, dtStep, parE);
        stressL = recover_nodal_stress_axisym_viscoelastic(out.meshL, st.uL, st.uLPrev, dtStep, parL);
    else
        stressE = recover_nodal_stress_axisym(out.meshE, st.uE, parE);
        stressL = recover_nodal_stress_axisym(out.meshL, st.uL, parL);
    end

    fluidColumns = [1, 3, 2, 4];
    componentTitles = {'\sigma_{rr}', '\sigma_{zz}', '\sigma_{\theta\theta}', '\sigma_{rz}'};

    solidStressE = {stressE.sigma_rr, stressE.sigma_zz, stressE.sigma_tt, stressE.sigma_rz};
    solidStressL = {stressL.sigma_rr, stressL.sigma_zz, stressL.sigma_tt, stressL.sigma_rz};

    parentObj = [];
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            if isgraphics(varargin{i}, 'figure') || isgraphics(varargin{i}, 'tiledlayout')
                parentObj = varargin{i};
                break;
            end
        end
    end

    if isempty(parentObj)
        fig = figure('Visible', 'off');
        tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    elseif isgraphics(parentObj, 'figure')
        tl = tiledlayout(parentObj, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    else
        tl = parentObj;
    end

    for k = 1:4
        ax = nexttile(tl, k);
        set(ax, 'FontSize', 16);
        hold(ax, 'on');

        selectedFluidStress = fluidSigma(:, fluidColumns(k));
        render_fluid_stress_component(ax, out, plotstep, fluidCenter, selectedFluidStress);

        plot_select_nodal_stress_contour(out.meshE, st.uE, solidStressE{k});
        plot_select_nodal_stress_contour(out.meshL, st.uL, solidStressL{k});

        draw_select_fluid_boundaries(ax, out, plotstep);
        colorbar(ax);

        xlabel(ax, 'r [\mum]');
        ylabel(ax, 'z [\mum]');
        title(ax, sprintf('%s, t = %.4g s', componentTitles{k}, out.t(plotstep)));

        axis(ax, 'equal', 'tight');
        rMaxE = max(out.meshE.nodes(:,1)) * 1e6;
        xlim(ax, [0 rMaxE]);
        box(ax, 'on');
    end

    if useViscoRecovery
        title(tl, 'Solid (elastic + Kelvin-Voigt viscous) and fluid (total) stress components');
    else
        title(tl, 'Solid (elastic ONLY) and fluid (total) stress components');
    end
end

function render_fluid_stress_component(ax, out, plotstep, fluidCenter, fluidField)
    if isfield(out, 'native2D') && isfield(out.native2D, 'R') && ...
            size(out.native2D.R, 1) > 1 && size(out.native2D.R, 2) > 1
        R = out.native2D.R * 1e6;
        Z = out.native2D.Z * 1e6;
        if numel(fluidField) == numel(R)
            S = reshape(fluidField, size(R));
            contourf(ax, R, Z, S, 48, 'LineColor', 'none');
            return;
        end
    end

    native2D.S = fluidField(:);
    native2D.R = fluidCenter(:,1);
    native2D.Z = fluidCenter(:,2);

    valid = isfinite(native2D.R) & isfinite(native2D.Z) & isfinite(native2D.S);
    if nnz(valid) < 3, return; end

    rVals = native2D.R(valid);
    zVals = native2D.Z(valid);
    fVals = native2D.S(valid);

    rGrid = linspace(min(rVals), max(rVals), 200);
    zGrid = linspace(min(zVals), max(zVals), 200);
    [RGrid, ZGrid] = meshgrid(rGrid, zGrid);

    F = scatteredInterpolant(rVals(:), zVals(:), fVals(:), 'linear', 'none');
    FGrid = F(RGrid, ZGrid);

    statePlot = state_for_plot_at_step(out, plotstep);
    solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid);
    FGrid(solidMask) = NaN;

    finiteVals = FGrid(isfinite(FGrid));
    if isempty(finiteVals), return; end

    if max(finiteVals) > min(finiteVals)
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 48, 'LineColor', 'none');
    else
        contourf(ax, RGrid*1e6, ZGrid*1e6, FGrid, 1, 'LineColor', 'none');
    end
end

function solidMask = deformed_solids_mask_on_grid(out, statePlot, RGrid, ZGrid)
    solidMask = false(size(RGrid));
    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL)
        solidMask = solidMask | local_deformed_solid_mask(out.meshL, statePlot.uL, RGrid, ZGrid);
    end
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        solidMask = solidMask | local_deformed_solid_mask(out.meshE, statePlot.uE, RGrid, ZGrid);
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
        inBox = RGrid >= min(rv) & RGrid <= max(rv) & ZGrid >= min(zv) & ZGrid <= max(zv);
        if any(inBox(:))
            localMask = false(size(solidMask));
            localMask(inBox) = inpolygon(RGrid(inBox), ZGrid(inBox), rv, zv);
            solidMask = solidMask | localMask;
        end
    end
end
