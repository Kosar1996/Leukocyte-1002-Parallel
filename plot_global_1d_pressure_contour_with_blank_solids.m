function plot_global_1d_pressure_contour_with_blank_solids(out, par)
    if ~isfield(out, 'global1D') || ~isfield(out.global1D, 'PfullFinal')
        return;
    end

    global1D = out.global1D;
    r = global1D.r(:);
    z = global1D.z(:);
    if isfield(global1D, 'PblankFinal') && ~isempty(global1D.PblankFinal)
        Pplot = global1D.PblankFinal;
    else
        Pplot = global1D.PfullFinal;
    end
    if isempty(Pplot) || ~isequal(size(Pplot), [numel(r), numel(z)])
        warning('Skipping final global pressure contour: pressure array size does not match r-z grid.');
        return;
    end

    pVals = Pplot(isfinite(Pplot));
    if isempty(pVals)
        warning('Skipping final global pressure contour: no finite fluid pressure values remain after solid masking.');
        return;
    end

    figure;
    ax = gca;
    set(ax, 'FontSize', 24);
    if max(pVals) > min(pVals)
        contourf(ax, r*1e6, z*1e6, Pplot.', linspace(min(pVals), max(pVals), 41), ...
            'LineColor', 'none');
    else
        contourf(ax, r*1e6, z*1e6, Pplot.', 1, 'LineColor', 'none');
    end
    colorbar;
    hold(ax, 'on');

    rMin = min(r);
    rMax = max(r);
    zMin = min(z);
    zMax = max(z);
    plot(ax, [rMin rMax]*1e6, [zMin zMin]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [rMin rMax]*1e6, [zMax zMax]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [rMax rMax]*1e6, [zMin zMax]*1e6, 'Color', [0.10 0.45 0.95], 'LineWidth', 1.2);
    plot(ax, [rMin rMin]*1e6, [zMin zMax]*1e6, '--', 'Color', [0.10 0.45 0.95], 'LineWidth', 1.0);

    statePlot = out.state;
    if isfield(out, 'stateHist') && isfield(out, 'stopStep') && ...
            out.stopStep >= 1 && numel(out.stateHist) >= out.stopStep && ...
            ~isempty(out.stateHist{out.stopStep})
        statePlot = out.stateHist{out.stopStep};
    end

    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL) && ...
            ~(isfield(par, 'noLeukocyte') && par.noLeukocyte)
        draw_deformed_solid_blank(out.meshL, statePlot.uL, [0.45 0.65 0.45]);
        if isfield(out, 'interfaceL') && ~isempty(out.interfaceL)
            [rLDef, zLDef] = deformed_interface_curve(out.meshL, statePlot.uL, out.interfaceL);
            plot(ax, rLDef*1e6, zLDef*1e6, 'Color', [0.00 0.35 0.10], 'LineWidth', 1.8);
        end
    end
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        draw_deformed_solid_blank(out.meshE, statePlot.uE, [0.55 0.62 0.72]);
        if isfield(out, 'interfaceE') && ~isempty(out.interfaceE)
            [rEDef, zEDef] = deformed_interface_curve(out.meshE, statePlot.uE, out.interfaceE);
            plot(ax, rEDef*1e6, zEDef*1e6, 'Color', [0.00 0.15 0.65], 'LineWidth', 1.8);
        end
    end

    axis(ax, 'equal');
    xlim(ax, [rMin rMax]*1e6);
    ylim(ax, [zMin zMax]*1e6);
    box(ax, 'on');
    grid(ax, 'off');
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Final global pressure contour, solids blank, t = %.4g s', out.t(out.stopStep)));
end

function h = draw_deformed_solid_blank(mesh, u, edgeColor)
    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    h = patch('Faces', mesh.conn, ...
        'Vertices', [rDef, zDef]*1e6, ...
        'FaceColor', 'w', ...
        'EdgeColor', edgeColor, ...
        'LineWidth', 0.15);
end