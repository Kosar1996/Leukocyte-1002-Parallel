%% =========================================================================
% AUDITED NATIVE 2D VELOCITY PLOTTING FUNCTION (OPTION A: AUTO-FILL)
% File & Function Name: plot_select_native2d_velocity.m
% Fills full 2D velocity domain using nearest-neighbor far-field extrapolation.
% =========================================================================

function plot_select_native2d_velocity(out, plotstep, varargin)
    if nargin < 2 || isempty(plotstep)
        if isfield(out, 'stopStep') && ~isempty(out.stopStep)
            plotstep = out.stopStep;
        else
            plotstep = 1;
        end
    end

    ax = [];
    if ~isempty(varargin)
        for i = 1:numel(varargin)
            if isgraphics(varargin{i}, 'axes')
                ax = varargin{i};
                break;
            elseif isgraphics(varargin{i}, 'figure')
                ax = gca(varargin{i});
                break;
            end
        end
    end
    if isempty(ax)
        fig = figure('Visible', 'off');
        ax = gca(fig);
    end

    set(ax, 'FontSize', 22);

    if isfield(out, 't') && ~isempty(out.t)
        if numel(out.t) >= plotstep && plotstep > 0
            tVal = out.t(plotstep);
        else
            tVal = out.t(end);
        end
    else
        tVal = 0.0;
    end

    [R_raw, Z_raw, Uz_raw] = extract_velocity_fields(out, plotstep);

    if isempty(Uz_raw) || ~any(isfinite(Uz_raw(:)))
        warning('plot_select_native2d_velocity: No finite velocity field available.');
        return;
    end

    valid = isfinite(R_raw) & isfinite(Z_raw) & isfinite(Uz_raw);
    rVals = R_raw(valid) * 1e6;
    zVals = Z_raw(valid) * 1e6;
    uzVals = Uz_raw(valid) * 1e6; % um/s

    rMin = 0; rMax = max(rVals); if rMax <= 0, rMax = 6; end
    zMin = min(zVals); zMax = max(zVals);
    if zMin == zMax, zMin = -6; zMax = 10; end

    [RGrid, ZGrid] = meshgrid(linspace(rMin, rMax, 250), linspace(zMin, zMax, 250));

    F = scatteredInterpolant(rVals(:), zVals(:), uzVals(:), 'linear', 'nearest');
    UzGrid = F(RGrid, ZGrid);

    contourf(ax, RGrid, ZGrid, UzGrid, 48, 'LineColor', 'none');
    hold(ax, 'on');

    statePlot = struct();
    if isfield(out, 'stateHist') && numel(out.stateHist) >= plotstep && ~isempty(out.stateHist{plotstep})
        statePlot = out.stateHist{plotstep};
    elseif isfield(out, 'state')
        statePlot = out.state;
    end
    draw_filled_solid_cells(ax, out, statePlot);

    if exist('draw_select_solid_blanks', 'file') == 2
        try draw_select_solid_blanks(ax, out, plotstep); catch; end
    end
    if exist('draw_select_fluid_boundaries', 'file') == 2
        try draw_select_fluid_boundaries(ax, out, plotstep); catch; end
    end

    colorbar(ax);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Axial velocity u_z(r,z), t = %.4g s', tVal));
    xlim(ax, [rMin rMax]);
    ylim(ax, [zMin zMax]);
    axis(ax, 'equal');
    box(ax, 'on');
end

function [R, Z, Uz] = extract_velocity_fields(out, plotstep)
    R = []; Z = []; Uz = [];
    if isfield(out, 'native2D') && isstruct(out.native2D) && isfield(out.native2D, 'uz')
        Uz = out.native2D.uz;
        if isfield(out.native2D, 'R'), R = out.native2D.R; end
        if isfield(out.native2D, 'Z'), Z = out.native2D.Z; end
    end

    if isempty(Uz) && isfield(out, 'uzCHist') && ~isempty(out.uzCHist)
        if ndims(out.uzCHist) == 3 && size(out.uzCHist, 3) >= plotstep
            Uz = out.uzCHist(:,:,plotstep);
            if isfield(out, 'RPHist') && size(out.RPHist,3) >= plotstep, R = out.RPHist(:,:,plotstep); end
            if isfield(out, 'ZPHist') && size(out.ZPHist,3) >= plotstep, Z = out.ZPHist(:,:,plotstep); end
        end
    end

    if ~isempty(Uz) && (isempty(R) || isempty(Z) || ~isequal(size(R), size(Uz)))
        [nR, nZ] = size(Uz);
        [R, Z] = meshgrid(linspace(0, 6e-6, nZ), linspace(-6e-6, 10e-6, nR));
    end
end

function draw_filled_solid_cells(ax, out, statePlot)
    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ~isempty(out.meshL) && ~isempty(statePlot.uL)
        rDefL = (out.meshL.nodes(:,1) + statePlot.uL(1:2:end)) * 1e6;
        zDefL = (out.meshL.nodes(:,2) + statePlot.uL(2:2:end)) * 1e6;
        if isfield(out.meshL, 'conn')
            for e = 1:size(out.meshL.conn, 1)
                ids = out.meshL.conn(e,:);
                patch(ax, rDefL(ids), zDefL(ids), [0.85 0.85 0.85], 'EdgeColor', [0.2 0.2 0.2]);
            end
        end
    end
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ~isempty(out.meshE) && ~isempty(statePlot.uE)
        rDefE = (out.meshE.nodes(:,1) + statePlot.uE(1:2:end)) * 1e6;
        zDefE = (out.meshE.nodes(:,2) + statePlot.uE(2:2:end)) * 1e6;
        if isfield(out.meshE, 'conn')
            for e = 1:size(out.meshE.conn, 1)
                ids = out.meshE.conn(e,:);
                patch(ax, rDefE(ids), zDefE(ids), [0.75 0.75 0.75], 'EdgeColor', [0.2 0.2 0.2]);
            end
        end
    end
end
