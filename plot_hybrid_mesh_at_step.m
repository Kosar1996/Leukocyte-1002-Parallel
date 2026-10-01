function plot_hybrid_mesh_at_step(out, hybridMesh, stepIdx)
    if nargin < 3 || isempty(stepIdx)
        stepIdx = 1;
    end
    if isempty(fieldnames(hybridMesh))
        warning('No hybrid mesh is available at step %d to plot.', stepIdx);
        return;
    end

    statePlot = state_for_plot_at_step(out, stepIdx);
    fluidPlot = [];
    if isfield(out, 'fluidHist') && numel(out.fluidHist) >= stepIdx
        fluidPlot = out.fluidHist{stepIdx};
    end

    figure('Name', sprintf('Hybrid mesh step %d', stepIdx), 'Color', 'w');
    ax = gca;
    hold(ax, 'on');
    set(ax, 'FontSize', 18);

    hLegend = gobjects(0);
    [drewRuntimeBlocks, hRuntime] = draw_runtime_hybrid_blocks(ax, fluidPlot);
    hLegend = [hLegend; hRuntime(:)];
    if ~drewRuntimeBlocks
        hUp = draw_hybrid_summary_block( ...
            ax, hybridMesh.exterior2D.upstream, ...
            [0.40 0.55 0.75], 0.45, '2D exterior mesh');
        hDown = draw_hybrid_summary_block( ...
            ax, hybridMesh.exterior2D.downstream, ...
            [0.40 0.55 0.75], 0.45, '2D exterior mesh');
        if ~isempty(hUp) && isgraphics(hUp)
            hLegend = [hLegend; hUp];
            if ~isempty(hDown) && isgraphics(hDown)
                set(hDown, 'HandleVisibility', 'off');
            end
        elseif ~isempty(hDown) && isgraphics(hDown)
            hLegend = [hLegend; hDown];
        end
    end

    hGap = draw_hybrid_gap_mesh(ax, hybridMesh);
    hLegend = [hLegend; hGap(:)];
    hBoundary = draw_fluid_boundaries_for_state(ax, out, statePlot);
    hLegend = [hLegend; hBoundary(:)];

    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    if isfield(hybridMesh, 't') && isfinite(hybridMesh.t)
        title(ax, sprintf('Hybrid mesh at first time step: step %d, t = %.4g s', ...
            stepIdx, hybridMesh.t));
    else
        title(ax, sprintf('Hybrid mesh at first time step: step %d', stepIdx));
    end
    axis(ax, 'equal');
    axis(ax, 'tight');
    box(ax, 'on');

    hLegend = hLegend(isgraphics(hLegend));
    if ~isempty(hLegend)
        legend(ax, hLegend, get(hLegend, 'DisplayName'), ...
            'Location', 'best');
    end
end

function [drew, hLegend] = draw_runtime_hybrid_blocks(ax, fluidPlot)
    drew = false;
    hLegend = gobjects(0);
    if ~isstruct(fluidPlot) || ~isfield(fluidPlot, 'meshF') || ...
            ~isfield(fluidPlot.meshF, 'blocks')
        return;
    end

    blocks = fluidPlot.meshF.blocks;
    hUp = draw_runtime_hybrid_block(ax, blocks.upstream, ...
        [0.40 0.55 0.75], 0.45, '2D exterior mesh');
    hDown = draw_runtime_hybrid_block(ax, blocks.downstream, ...
        [0.40 0.55 0.75], 0.45, '2D exterior mesh');
    hCenters = draw_runtime_pressure_centers(ax, blocks);

    drew = (~isempty(hUp) && isgraphics(hUp)) || ...
           (~isempty(hDown) && isgraphics(hDown));
    if ~isempty(hUp) && isgraphics(hUp)
        hBlock = hUp;
        if ~isempty(hDown) && isgraphics(hDown)
            set(hDown, 'HandleVisibility', 'off');
        end
    else
        hBlock = hDown;
    end

    hNew = [hBlock; hCenters];
    hNew = hNew(isgraphics(hNew));
    if ~isempty(hNew)
        hLegend = hNew(:);
    end
end

function hFirst = draw_runtime_hybrid_block(ax, block, color, lineWidth, displayName)
    hFirst = gobjects(0);
    if ~isstruct(block) || ~isfield(block, 'mesh') || ~isstruct(block.mesh)
        return;
    end

    mesh = block.mesh;
    if isfield(mesh, 'Rzf') && isfield(mesh, 'zF') && ~isempty(mesh.Rzf)
        R = mesh.Rzf;
        Z = repmat(mesh.zF(:).', size(R, 1), 1);
    elseif isfield(mesh, 'Rp') && isfield(mesh, 'Zp') && ~isempty(mesh.Rp)
        R = mesh.Rp;
        Z = mesh.Zp;
    else
        return;
    end

    hLines = plot_mesh_matrix(ax, R, Z, color, lineWidth);
    if isempty(hLines)
        return;
    end
    hFirst = hLines(1);
    set(hFirst, 'DisplayName', displayName);
    if numel(hLines) > 1
        set(hLines(2:end), 'HandleVisibility', 'off');
    end
end

function hFirst = draw_runtime_pressure_centers(ax, blocks)
    hFirst = gobjects(0);
    blockNames = {'upstream', 'downstream'};
    for i = 1:numel(blockNames)
        name = blockNames{i};
        if ~isfield(blocks, name) || ~isstruct(blocks.(name)) || ...
                ~isfield(blocks.(name), 'mesh')
            continue;
        end
        mesh = blocks.(name).mesh;
        if ~isfield(mesh, 'Rp') || ~isfield(mesh, 'Zp') || isempty(mesh.Rp)
            continue;
        end
        h = plot(ax, mesh.Rp(:)*1e6, mesh.Zp(:)*1e6, '.', ...
            'Color', [0.10 0.25 0.80], 'MarkerSize', 4, ...
            'DisplayName', 'pressure centers');
        if isempty(hFirst)
            hFirst = h;
        else
            set(h, 'HandleVisibility', 'off');
        end
    end
end

function h = draw_hybrid_gap_mesh(ax, hybridMesh)
    h = gobjects(0);
    if ~isfield(hybridMesh, 'gap1D') || ~isfield(hybridMesh.gap1D, 'z') || ...
            isempty(hybridMesh.gap1D.z)
        return;
    end

    gap = hybridMesh.gap1D;
    h(end+1,1) = plot(ax, gap.rL*1e6, gap.z*1e6, '-', ...
        'Color', [0.05 0.05 0.05], 'LineWidth', 1.2, ...
        'DisplayName', '1D gap boundaries');
    hEnd = plot(ax, gap.rE*1e6, gap.z*1e6, '-', ...
        'Color', [0.05 0.05 0.05], 'LineWidth', 1.2);
    set(hEnd, 'HandleVisibility', 'off');

    h(end+1,1) = plot(ax, gap.rMid*1e6, gap.z*1e6, 'o--', ...
        'Color', [0.85 0.20 0.10], 'MarkerFaceColor', [0.85 0.20 0.10], ...
        'MarkerSize', 3.5, 'LineWidth', 0.9, ...
        'DisplayName', '1D gap nodes');

    if isfield(hybridMesh, 'gapZ') && isfield(hybridMesh, 'rDomain')
        rSpan = hybridMesh.rDomain(:).';
        for i = 1:numel(hybridMesh.gapZ)
            hCouple = plot(ax, rSpan*1e6, ...
                hybridMesh.gapZ(i)*[1 1]*1e6, ':', ...
                'Color', [0.15 0.15 0.15], 'LineWidth', 1.0);
            set(hCouple, 'HandleVisibility', 'off');
        end
    end
end

function h = draw_fluid_boundaries_for_state(ax, out, statePlot)
    h = gobjects(0);
    if ~isfield(out, 'z') || isempty(out.z)
        return;
    end
    z = out.z(:);

    if isfield(out, 'meshL') && isfield(out, 'interfaceL') && ...
            isfield(statePlot, 'uL') && ~isempty(statePlot.uL)
        [rLDef, zLDef] = local_deformed_interface_curve( ...
            out.meshL, statePlot.uL, out.interfaceL);
        h(end+1,1) = plot(ax, rLDef*1e6, zLDef*1e6, 'k-', ...
            'LineWidth', 1.4, 'DisplayName', 'solid/fluid boundaries');
    elseif isfield(statePlot, 'deltaL')
        h(end+1,1) = plot(ax, statePlot.deltaL(:)*1e6, z*1e6, 'k-', ...
            'LineWidth', 1.4, 'DisplayName', 'solid/fluid boundaries');
    end

    if isfield(out, 'meshE') && isfield(out, 'interfaceE') && ...
            isfield(statePlot, 'uE') && ~isempty(statePlot.uE)
        [rEDef, zEDef] = local_deformed_interface_curve( ...
            out.meshE, statePlot.uE, out.interfaceE);
        hE = plot(ax, rEDef*1e6, zEDef*1e6, 'k-', 'LineWidth', 1.4);
        set(hE, 'HandleVisibility', 'off');
    elseif isfield(statePlot, 'deltaE')
        hE = plot(ax, statePlot.deltaE(:)*1e6, z*1e6, 'k-', ...
            'LineWidth', 1.4);
        set(hE, 'HandleVisibility', 'off');
    end
end

function hFirst = draw_hybrid_summary_block(ax, block, color, lineWidth, displayName)
    hFirst = gobjects(0);
    if ~isstruct(block) || ~isfield(block, 'R') || ~isfield(block, 'Z') || ...
            isempty(block.R)
        return;
    end

    hLines = plot_mesh_matrix(ax, block.R, block.Z, color, lineWidth);
    if isempty(hLines)
        return;
    end
    hFirst = hLines(1);
    set(hFirst, 'DisplayName', displayName);
    if numel(hLines) > 1
        set(hLines(2:end), 'HandleVisibility', 'off');
    end
end



