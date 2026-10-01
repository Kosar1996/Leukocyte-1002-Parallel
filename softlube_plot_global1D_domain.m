function diag = softlube_plot_global1D_domain(cfg, out, opts)
%SOFTLUBE_PLOT_GLOBAL1D_DOMAIN Plot the global 1D pressure domain and mesh.
%
% Usage:
%   cfg = softlube_case_global1D_deformable_leukocyte_P600();
%   diag = softlube_plot_global1D_domain(cfg);
%
%   out = softlube_run_case_global1D_coupled(cfg);
%   diag = softlube_plot_global1D_domain(cfg, out);

    if nargin < 1 || isempty(cfg)
        cfg = softlube_case_global1D_deformable_leukocyte_P600();
    end
    if nargin < 2
        out = [];
    end
    if nargin < 3 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'makeFigure')
        opts.makeFigure = true;
    end
    if ~isfield(opts, 'showTractionSupport')
        opts.showTractionSupport = false;
    end

    z = global1d_z_grid_from_cfg(cfg);
    rDomain = [0, 15e-6];
    if isfield(cfg, 'global1D') && isfield(cfg.global1D, 'rDomain')
        rDomain = cfg.global1D.rDomain;
    end

    SE = load(cfg.geometry.endotheliumPrestressFile);
    SL = load(cfg.geometry.leukocytePrestressFile);

    uE = SE.uE_pre(:);
    uL = SL.uL_pre(:);
    if ~isempty(out) && isstruct(out) && isfield(out, 'state')
        if isfield(out.state, 'uE') && ~isempty(out.state.uE)
            uE = out.state.uE(:);
        end
        if isfield(out.state, 'uL') && ~isempty(out.state.uL)
            uL = out.state.uL(:);
        end
    end

    [rE, zE] = interface_curve(SE.meshE, uE, SE.interfaceE);
    [rL, zL] = interface_curve(SL.meshL, uL, SL.interfaceL);
    [VE, FE] = deformed_mesh_patch(SE.meshE, uE);
    [VL, FL] = deformed_mesh_patch(SL.meshL, uL);

    rBaseE = [];
    zBaseE = [];
    if isfield(SE, 'baseE')
        [rBaseE, zBaseE] = interface_curve(SE.meshE, uE, SE.baseE);
    end
    rBaseL = [];
    zBaseL = [];
    if isfield(SL, 'baseL')
        [rBaseL, zBaseL] = interface_curve(SL.meshL, uL, SL.baseL);
    end

    fineWindow = [-2e-6, 6e-6];
    if isfield(cfg, 'parOverrides') && isfield(cfg.parOverrides, 'global1DFineWindow')
        fineWindow = cfg.parOverrides.global1DFineWindow;
    end

    diag = struct();
    diag.z = z;
    diag.dz = diff(z);
    diag.rDomain = rDomain;
    diag.fineWindow = fineWindow;
    diag.endothelium = struct('r', rE, 'z', zE);
    diag.leukocyte = struct('r', rL, 'z', zL);
    diag.endotheliumOuterBoundary = struct('r', rBaseE, 'z', zBaseE);
    diag.leukocyteAxisBoundary = struct('r', rBaseL, 'z', zBaseL);
    diag.numNodes = numel(z);
    diag.minDz = min(diff(z));
    diag.maxDz = max(diff(z));
    diag.pressureBC = struct( ...
        'zMin', cfg.fluid.zMin, ...
        'zMax', cfg.fluid.zMax, ...
        'pIn', cfg.bc.pIn, ...
        'pOut', cfg.bc.pOut, ...
        'rMaxPressure', 0, ...
        'axisCondition', 'dP/dr = 0');
    if isfield(cfg, 'interface') && isfield(cfg.interface, 'leukocytePressureSupportZ')
        diag.leukocytePressureSupportZ = cfg.interface.leukocytePressureSupportZ(:).';
    end

    if ~opts.makeFigure
        return;
    end

    figure;
    hold on;
    set(gca, 'FontSize', 18);

    zMin = cfg.fluid.zMin;
    zMax = cfg.fluid.zMax;
    rMin = rDomain(1);
    rMax = rDomain(2);

    hFine = patch([rMin, rMax, rMax, rMin]*1e6, ...
          [fineWindow(1), fineWindow(1), fineWindow(2), fineWindow(2)]*1e6, ...
          [0.92 0.96 1.00], ...
          'EdgeColor', 'none', 'FaceAlpha', 0.75);

    rectangle('Position', [rMin*1e6, zMin*1e6, ...
        (rMax-rMin)*1e6, (zMax-zMin)*1e6], ...
        'EdgeColor', [0.1 0.1 0.1], 'LineWidth', 1.5, ...
        'HandleVisibility', 'off');
    hBox = plot(nan, nan, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);

    hMesh = plot(nan, nan, '-', 'Color', [0.55 0.70 0.90], 'LineWidth', 0.7);
    for k = 1:numel(z)
        if z(k) < fineWindow(1) || z(k) > fineWindow(2)
            color = [0.75 0.75 0.75];
            lw = 0.5;
        else
            color = [0.55 0.70 0.90];
            lw = 0.7;
        end
        plot([rMin, rMax]*1e6, [z(k), z(k)]*1e6, '-', ...
            'Color', color, 'LineWidth', lw, 'HandleVisibility', 'off');
    end

    hLeuSolid = patch('Faces', FL, 'Vertices', VL*1e6, ...
        'FaceColor', [0.70 0.84 1.00], 'FaceAlpha', 0.35, ...
        'EdgeColor', [0.60 0.75 0.95], 'LineWidth', 0.2);
    hEndoSolid = patch('Faces', FE, 'Vertices', VE*1e6, ...
        'FaceColor', [1.00 0.78 0.72], 'FaceAlpha', 0.35, ...
        'EdgeColor', [0.95 0.65 0.60], 'LineWidth', 0.2);

    hLeuInterface = plot(rL*1e6, zL*1e6, ...
        'Color', [0.05 0.28 0.75], 'LineWidth', 2.4);
    hEndoInterface = plot(rE*1e6, zE*1e6, ...
        'Color', [0.78 0.12 0.08], 'LineWidth', 2.4);

    hPz = plot([rMin rMax]*1e6, [zMin zMin]*1e6, '--', ...
        'Color', [0.65 0 0], 'LineWidth', 2);
    plot([rMin rMax]*1e6, [zMax zMax]*1e6, '--', ...
        'Color', [0.65 0 0], 'LineWidth', 2, 'HandleVisibility', 'off');

    endoOuterBlock = [];
    if ~isempty(zBaseE)
        endoOuterBlock = [min(zBaseE), max(zBaseE)];
    end
    hPrMax = plot_boundary_segments(rMax, [zMin zMax], endoOuterBlock, ...
        '-', [0.65 0 0], 2);

    leuAxisBlock = [];
    if ~isempty(zBaseL)
        leuAxisBlock = [min(zBaseL), max(zBaseL)];
    end
    hAxisFluid = plot_boundary_segments(0, [zMin zMax], leuAxisBlock, ...
        '-', [0.10 0.30 0.70], 2);

    hLeuAxis = [];
    if ~isempty(rBaseL)
        hLeuAxis = plot(rBaseL*1e6, zBaseL*1e6, ...
            'Color', [0.00 0.15 0.50], 'LineWidth', 2.6);
    end
    hEndoOuter = [];
    if ~isempty(rBaseE)
        hEndoOuter = plot(rBaseE*1e6, zBaseE*1e6, ...
            'Color', [0.50 0.00 0.00], 'LineWidth', 2.6);
    end

    showSupportLegend = false;
    hSupport = [];
    if opts.showTractionSupport && isfield(cfg, 'interface') && ...
            isfield(cfg.interface, 'leukocytePressureSupportZ')
        zs = cfg.interface.leukocytePressureSupportZ(:).';
        supportR = rMax + 0.35e-6;
        hSupport = plot([supportR supportR]*1e6, zs*1e6, '--', ...
            'Color', [0.15 0.55 0.25], 'LineWidth', 2);
        showSupportLegend = true;
    end

    xlabel('r [\mum]');
    ylabel('z [\mum]');
    title(sprintf('Global 1D pressure domain: Nz=%d, dz=[%.2f, %.2f] \\mum', ...
        numel(z), min(diff(z))*1e6, max(diff(z))*1e6));
    legendHandles = [ ...
        plot(nan, nan, 's', 'MarkerFaceColor', [0.92 0.96 1.00], ...
            'MarkerEdgeColor', [0.92 0.96 1.00], 'MarkerSize', 12), ...
        plot(nan, nan, '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5), ...
        plot(nan, nan, '-', 'Color', [0.55 0.70 0.90], 'LineWidth', 0.7), ...
        plot(nan, nan, 's', 'MarkerFaceColor', [0.70 0.84 1.00], ...
            'MarkerEdgeColor', [0.60 0.75 0.95], 'MarkerSize', 12), ...
        plot(nan, nan, 's', 'MarkerFaceColor', [1.00 0.78 0.72], ...
            'MarkerEdgeColor', [0.95 0.65 0.60], 'MarkerSize', 12), ...
        plot(nan, nan, '-', 'Color', [0.05 0.28 0.75], 'LineWidth', 2.4), ...
        plot(nan, nan, '-', 'Color', [0.78 0.12 0.08], 'LineWidth', 2.4), ...
        plot(nan, nan, '--', 'Color', [0.65 0 0], 'LineWidth', 2), ...
        plot(nan, nan, '-', 'Color', [0.65 0 0], 'LineWidth', 2), ...
        plot(nan, nan, '-', 'Color', [0.10 0.30 0.70], 'LineWidth', 2)];
    legendEntries = {'fine pressure window', 'global box', 'axial pressure mesh', ...
        'leukocyte solid', 'endothelium solid', ...
        'leukocyte fluid interface', 'endothelium fluid interface', ...
        'P=0 at z=-6/10', 'P=0 at exposed r=15', ...
        'dP/dr=0 at exposed r=0'};
    if ~isempty(hLeuAxis)
        legendHandles(end+1) = plot(nan, nan, '-', ...
            'Color', [0.00 0.15 0.50], 'LineWidth', 2.6); %#ok<AGROW>
        legendEntries{end+1} = 'leukocyte axis solid support'; %#ok<AGROW>
    end
    if ~isempty(hEndoOuter)
        legendHandles(end+1) = plot(nan, nan, '-', ...
            'Color', [0.50 0.00 0.00], 'LineWidth', 2.6); %#ok<AGROW>
        legendEntries{end+1} = 'endothelium outer solid boundary'; %#ok<AGROW>
    end
    if showSupportLegend
        legendHandles(end+1) = plot(nan, nan, '--', ...
            'Color', [0.15 0.55 0.25], 'LineWidth', 2); %#ok<AGROW>
        legendEntries{end+1} = 'leukocyte traction support';
    end
    legend(legendHandles, legendEntries, 'Location', 'eastoutside');
    if showSupportLegend
        xlim([rMin, rMax + 0.7e-6]*1e6);
    else
        xlim([rMin, rMax]*1e6);
    end
    ylim([zMin, zMax]*1e6);
    axis normal;
    grid on;
end

function z = global1d_z_grid_from_cfg(cfg)
    if ~isfield(cfg, 'parOverrides') || ...
            ~isfield(cfg.parOverrides, 'useGlobal1DCoarseEdgeMesh') || ...
            ~cfg.parOverrides.useGlobal1DCoarseEdgeMesh
        Nz = 121;
        if isfield(cfg, 'parOverrides') && isfield(cfg.parOverrides, 'global1DNzFluid')
            Nz = cfg.parOverrides.global1DNzFluid;
        end
        z = linspace(cfg.fluid.zMin, cfg.fluid.zMax, Nz).';
        return;
    end

    fineWindow = cfg.parOverrides.global1DFineWindow;
    coarseDz = cfg.parOverrides.global1DCoarseDz;
    fineDz = cfg.parOverrides.global1DFineDz;
    fineWindow(1) = max(fineWindow(1), cfg.fluid.zMin);
    fineWindow(2) = min(fineWindow(2), cfg.fluid.zMax);

    z = sort([ ...
        segment_grid(cfg.fluid.zMin, fineWindow(1), coarseDz)
        segment_grid(fineWindow(1), fineWindow(2), fineDz)
        segment_grid(fineWindow(2), cfg.fluid.zMax, coarseDz)
        cfg.fluid.zMin
        cfg.fluid.zMax
        fineWindow(:)
        0
        4e-6]);

    z = z(z >= cfg.fluid.zMin - 1e-15 & z <= cfg.fluid.zMax + 1e-15);
    z(abs(z) < 1e-15) = 0;
    z(abs(z - 4e-6) < 1e-15) = 4e-6;
    z = z([true; diff(z) > 1e-15]);
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

function [V, F] = deformed_mesh_patch(mesh, u)
    r = mesh.nodes(:,1) + u(1:2:end);
    z = mesh.nodes(:,2) + u(2:2:end);
    V = [r(:), z(:)];
    F = mesh.conn;
end

function h = plot_boundary_segments(r, zDomain, blocked, lineStyle, color, lineWidth)
    h = plot(nan, nan, lineStyle, 'Color', color, 'LineWidth', lineWidth);

    segments = interval_minus(zDomain(:).', blocked);
    for k = 1:size(segments,1)
        if segments(k,2) > segments(k,1)
            plot([r r]*1e6, segments(k,:)*1e6, lineStyle, ...
                'Color', color, 'LineWidth', lineWidth, ...
                'HandleVisibility', 'off');
        end
    end
end

function segments = interval_minus(domain, blocked)
    domain = sort(domain(:).');
    if isempty(blocked)
        segments = domain;
        return;
    end

    blocked = sort(blocked(:).');
    blocked(1) = max(blocked(1), domain(1));
    blocked(2) = min(blocked(2), domain(2));

    segments = zeros(0,2);
    if blocked(1) > domain(1)
        segments(end+1,:) = [domain(1), blocked(1)]; %#ok<AGROW>
    end
    if blocked(2) < domain(2)
        segments(end+1,:) = [blocked(2), domain(2)]; %#ok<AGROW>
    end
    if isempty(segments)
        segments = [nan nan];
    end
end

function [rDef, zDef] = interface_curve(mesh, u, interfaceNodes)
    ids = interfaceNodes(:);
    rDef = mesh.nodes(ids,1) + u(2*ids - 1);
    zDef = mesh.nodes(ids,2) + u(2*ids);

    [zDef, order] = sort(zDef);
    rDef = rDef(order);
    [zUnique, ~, group] = unique(zDef);
    if numel(zUnique) < numel(zDef)
        rDef = accumarray(group, rDef, [], @mean);
        zDef = zUnique;
    end
end
