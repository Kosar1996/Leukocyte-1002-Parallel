function generate_heatmap_videos(matFile, label, outDir, fps)
%GENERATE_HEATMAP_VIDEOS Heatmap-vs-time videos (pressure/stress/velocity)
% for one case, replacing the final-timestep-only static images used in
% the Sep 9 report. Reuses the same plot_select_native2d_* helpers as
% generate_heatmaps.m, looping over every accepted step instead of just
% the last one.
%
% Usage:
%   generate_heatmap_videos('/path/to/out_case_x.mat', 'case_x', '/path/to/outdir', 5)
%
% matFile : path to a saved case output (.mat with variable 'out', or
%           'outRestart' for restart-style files)
% label   : case label used as the video filename prefix
% outDir  : output directory for the .mp4 files (created if missing)
% fps     : frames per second (default 5)

if nargin < 4 || isempty(fps)
    fps = 5;
end
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

set(0, 'DefaultFigureVisible', 'off');

S = load(matFile);
if isfield(S, 'out')
    out = S.out;
elseif isfield(S, 'outRestart')
    out = S.outRestart;
else
    error('generate_heatmap_videos:noOut', 'No ''out'' or ''outRestart'' variable found in %s', matFile);
end

nSteps = out.stopStep;
if nSteps < 1
    error('generate_heatmap_videos:noSteps', '%s has stopStep=%d, nothing to render', matFile, nSteps);
end

fields = {'pressure', @plot_select_native2d_pressure; ...
          'stress',   @plot_select_native2d_stress; ...
          'velocity', @plot_select_native2d_velocity};

for f = 1:size(fields,1)
    fieldName = fields{f,1};
    plotFn = fields{f,2};
    outFile = fullfile(outDir, [label '_' fieldName '.mp4']);
    fprintf('=== %s: %s ===\n', label, fieldName);

    try
        v = VideoWriter(outFile, 'MPEG-4');
        v.FrameRate = fps;
        open(v);

        framesWritten = 0;

        % t=0 reference frame (Sep 11): shows the initial prestressed shape
        % before any time-stepping, if the solver saved it (out.t0State/
        % out.t0Fluid -- runs predating this feature won't have these, and
        % are skipped gracefully rather than erroring).
        if isfield(out, 't0State') && isfield(out, 't0Fluid') && ~isempty(out.t0Fluid)
            try
                out0 = out;
                out0.fluidHist = {out.t0Fluid};
                out0.stateHist = {out.t0State};
                % is_hybrid_fluid_plot.m (used by all 3 plot_select_native2d_*
                % functions) gates on out.fluidHist{out.stopStep}, not on the
                % plotstep actually being drawn -- so stopStep must match this
                % 1-element synthetic history, or that check silently fails
                % and pressure's plot function falls into a dead fallback
                % branch that references an undefined variable. Likewise the
                % *Hist 3D arrays (PHist/RPHist/ZPHist) are what
                % plot_select_native2d_pressure.m actually contours from when
                % that gate passes -- left un-overridden they'd still hold
                % the real step-1 data, silently mislabeled as t=0.
                out0.stopStep = 1;
                out0.t(1) = 0; % FIX: label this frame's title as t=0, not the real step-1 time
                out0.PHist  = out.t0Fluid.P;
                out0.RPHist = out.t0Fluid.meshF.Rp;
                out0.ZPHist = out.t0Fluid.meshF.Zp;
                out0.native2D = struct();
                out0.native2D.P     = out.t0Fluid.P;
                out0.native2D.R     = out.t0Fluid.meshF.Rp;
                out0.native2D.Z     = out.t0Fluid.meshF.Zp;
                out0.native2D.ur    = out.t0Fluid.urC;
                out0.native2D.uz    = out.t0Fluid.uzC;
                out0.native2D.speed = sqrt(out.t0Fluid.urC.^2 + out.t0Fluid.uzC.^2);

                plotFn(out0, 1);
                frame = getframe(gcf);
                writeVideo(v, frame);
                close(gcf);
                framesWritten = framesWritten + 1;
            catch ME0
                fprintf('  t=0 frame: FAILED -- %s (skipped)\n', ME0.message);
                if exist('gcf', 'var') && ~isempty(get(0,'CurrentFigure'))
                    close(gcf);
                end
            end
        end

        for k = 1:nSteps
            try
                outK = out;
                outK.native2D = struct();
                outK.native2D.P     = out.PHist(:,:,k);
                outK.native2D.R     = out.RPHist(:,:,k);
                outK.native2D.Z     = out.ZPHist(:,:,k);
                outK.native2D.ur    = out.urCHist(:,:,k);
                outK.native2D.uz    = out.uzCHist(:,:,k);
                outK.native2D.speed = out.speedCHist(:,:,k);

                plotFn(outK, k);
                frame = getframe(gcf);
                writeVideo(v, frame);
                close(gcf);
                framesWritten = framesWritten + 1;
            catch MEk
                fprintf('  step %d: FAILED -- %s (skipped)\n', k, MEk.message);
                if exist('gcf', 'var') && ~isempty(get(0,'CurrentFigure'))
                    close(gcf);
                end
            end
        end

        close(v);
        hasT0Frame = isfield(out, 't0State') && isfield(out, 't0Fluid') && ~isempty(out.t0Fluid);
        expectedFrames = nSteps + double(hasT0Frame);
        fprintf('  %s: OK, %d/%d frames written (%d accepted step%s%s) -> %s\n', fieldName, ...
            framesWritten, expectedFrames, nSteps, ternary_local(nSteps==1,'','s'), ...
            ternary_local(hasT0Frame, ' + 1 t=0 frame', ''), outFile);
    catch ME
        fprintf('  %s: FAILED -- %s\n', fieldName, ME.message);
    end
end

fprintf('\nDone: %s\n', label);
end
