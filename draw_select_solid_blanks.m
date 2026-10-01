function draw_select_solid_blanks(ax, out,plotstep)
    statePlot = out.state;
    if isfield(out, 'stateHist') && out.stopStep >= 1 && ...
            numel(out.stateHist) >= out.stopStep && ~isempty(out.stateHist{plotstep})
        statePlot = out.stateHist{plotstep};
    end

    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ...
            ~isempty(out.meshL) && ~isempty(statePlot.uL)
        draw_deformed_solid_patch(ax, out.meshL, statePlot.uL, [0.05 0.35 0.12]);
    end

    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ...
            ~isempty(out.meshE) && ~isempty(statePlot.uE)
        draw_deformed_solid_patch(ax, out.meshE, statePlot.uE, [0.00 0.15 0.65]);
    end
end

function draw_deformed_solid_patch(ax, mesh, u, edgeColor)
    rDef = mesh.nodes(:,1) + u(1:2:end);
    zDef = mesh.nodes(:,2) + u(2:2:end);
    patch(ax, 'Faces', mesh.conn, ...
        'Vertices', [rDef, zDef]*1e6, ...
        'FaceColor', 'w', ...
        'EdgeColor', edgeColor, ...
        'LineWidth', 0.12);
end