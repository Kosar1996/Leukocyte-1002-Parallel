function draw_select_fluid_boundaries(ax, out,plotstep)
    statePlot = state_for_plot_at_step(out, plotstep);
    z = out.z(:);

    drewL = false;
    if isfield(out, 'meshL') && isfield(out, 'interfaceL') && ...
            isfield(statePlot, 'uL') && ~isempty(statePlot.uL)
        [rLDef, zLDef] = local_deformed_interface_curve( ...
            out.meshL, statePlot.uL, out.interfaceL);
        plot(ax, rLDef*1e6, zLDef*1e6, 'k-', 'LineWidth', 1.2);
        drewL = true;
    end
    if ~drewL && isfield(statePlot, 'deltaL')
        plot(ax, statePlot.deltaL(:)*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
    end

    drewE = false;
    if isfield(out, 'meshE') && isfield(out, 'interfaceE') && ...
            isfield(statePlot, 'uE') && ~isempty(statePlot.uE)
        [rEDef, zEDef] = local_deformed_interface_curve( ...
            out.meshE, statePlot.uE, out.interfaceE);
        plot(ax, rEDef*1e6, zEDef*1e6, 'k-', 'LineWidth', 1.2);
        drewE = true;
    end
    if ~drewE && isfield(statePlot, 'deltaE')
        plot(ax, statePlot.deltaE(:)*1e6, z*1e6, 'k-', 'LineWidth', 1.2);
    end

    rMax = max(out.native2D.R(:), [], 'omitnan');
    zMinPlot = min(z);
    zMaxPlot = max(z);
    plot(ax, [0 rMax]*1e6, [zMinPlot zMinPlot]*1e6, 'k--', 'LineWidth', 1.0);
    plot(ax, [0 rMax]*1e6, [zMaxPlot zMaxPlot]*1e6, 'k--', 'LineWidth', 1.0);
end
