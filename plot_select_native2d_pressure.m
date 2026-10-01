%% =========================================================================
% AUDITED TWO-PHASE PRESSURE PLOTTING FUNCTION
% File & Function Name: plot_select_native2d_pressure.m
% Resolves fluid white voids via auto-fill and restores original solid phase
% pressure rendering via recover_nodal_stress_axisym.
% =========================================================================

function plot_select_native2d_pressure(out, plotstep, varargin)
    if nargin < 2 || isempty(plotstep)
        if isfield(out, 'stopStep') && ~isempty(out.stopStep)
            plotstep = out.stopStep;
        else
            plotstep = 1;
        end
    end

    % Safe extraction of target axes handle
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

    % Safe time extraction
    if isfield(out, 't') && ~isempty(out.t)
        if numel(out.t) >= plotstep && plotstep > 0
            tVal = out.t(plotstep);
        else
            tVal = out.t(end);
        end
    else
        tVal = 0.0;
    end

    % 1. FLUID PHASE: Extract and render background pressure with nearest fill
    [R_raw, Z_raw, P_raw] = extract_pressure_fields(out, plotstep);

    if ~isempty(P_raw) && any(isfinite(P_raw(:)))
        valid = isfinite(R_raw) & isfinite(Z_raw) & isfinite(P_raw);
        rVals = R_raw(valid) * 1e6; % um
        zVals = Z_raw(valid) * 1e6;
        pVals = P_raw(valid);

        rMin = 0; rMax = max(rVals); if rMax <= 0, rMax = 6; end
        zMin = min(zVals); zMax = max(zVals);
        if zMin == zMax, zMin = -6; zMax = 10; end

        [RGrid, ZGrid] = meshgrid(linspace(rMin, rMax, 250), linspace(zMin, zMax, 250));
        F = scatteredInterpolant(rVals(:), zVals(:), pVals(:), 'linear', 'nearest');
        PGrid = F(RGrid, ZGrid);

        contourf(ax, RGrid, ZGrid, PGrid, 48, 'LineColor', 'none');
        hold(ax, 'on');
    end

    % 2. SOLID PHASE: Render hydrostatic pressure P = -(sig_rr + sig_zz + sig_tt)/3
    statePlot = struct();
    if isfield(out, 'stateHist') && numel(out.stateHist) >= plotstep && ~isempty(out.stateHist{plotstep})
        statePlot = out.stateHist{plotstep};
    elseif isfield(out, 'state')
        statePlot = out.state;
    end

    % Render Endothelium Solid Pressure
    if isfield(out, 'meshE') && isfield(statePlot, 'uE') && ~isempty(out.meshE) && ~isempty(statePlot.uE)
        try
            if exist('recover_nodal_stress_axisym', 'file') == 2
                stressE = recover_nodal_stress_axisym(out.meshE, statePlot.uE, out.par);
                pSolidE = -(stressE.sigma_rr + stressE.sigma_zz + stressE.sigma_tt) / 3;
                if exist('plot_select_nodal_stress_contour', 'file') == 2
                    plot_select_nodal_stress_contour(out.meshE, statePlot.uE, pSolidE);
                end
            end
        catch ME
            fprintf('  [Warning] Endothelium solid pressure rendering skipped: %s\n', ME.message);
        end
    end

    % Render Leukocyte Solid Pressure
    if isfield(out, 'meshL') && isfield(statePlot, 'uL') && ~isempty(out.meshL) && ~isempty(statePlot.uL)
        try
            if exist('recover_nodal_stress_axisym', 'file') == 2
                stressL = recover_nodal_stress_axisym(out.meshL, statePlot.uL, out.par);
                pSolidL = -(stressL.sigma_rr + stressL.sigma_zz + stressL.sigma_tt) / 3;
                if exist('plot_select_nodal_stress_contour', 'file') == 2
                    plot_select_nodal_stress_contour(out.meshL, statePlot.uL, pSolidL);
                end
            end
        catch ME
            fprintf('  [Warning] Leukocyte solid pressure rendering skipped: %s\n', ME.message);
        end
    end

    % 3. BOUNDARIES & OUTLINES
    if exist('draw_select_fluid_boundaries', 'file') == 2
        try draw_select_fluid_boundaries(ax, out, plotstep); catch; end
    end

    colorbar(ax);
    xlabel(ax, 'r [\mum]');
    ylabel(ax, 'z [\mum]');
    title(ax, sprintf('Fluid & Solid Pressure P(r,z), t = %.4g s', tVal));
    axis(ax, 'equal', 'tight');
    box(ax, 'on');
end

function [R, Z, P] = extract_pressure_fields(out, plotstep)
    R = []; Z = []; P = [];
    if isfield(out, 'native2D') && isstruct(out.native2D) && isfield(out.native2D, 'P')
        P = out.native2D.P;
        if isfield(out.native2D, 'R'), R = out.native2D.R; end
        if isfield(out.native2D, 'Z'), Z = out.native2D.Z; end
    end

    if isempty(P) && isfield(out, 'PHist') && ~isempty(out.PHist)
        if ndims(out.PHist) == 3 && size(out.PHist, 3) >= plotstep
            P = out.PHist(:,:,plotstep);
            if isfield(out, 'RPHist') && size(out.RPHist,3) >= plotstep, R = out.RPHist(:,:,plotstep); end
            if isfield(out, 'ZPHist') && size(out.ZPHist,3) >= plotstep, Z = out.ZPHist(:,:,plotstep); end
        end
    end

    if ~isempty(P) && (isempty(R) || isempty(Z) || ~isequal(size(R), size(P)))
        [nR, nZ] = size(P);
        [R, Z] = meshgrid(linspace(0, 6e-6, nZ), linspace(-6e-6, 10e-6, nR));
    end
end
