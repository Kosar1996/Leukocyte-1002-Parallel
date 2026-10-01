%% =========================================================================
% AUDITED & UNIFIED PARALLEL HEATMAP RENDERER (OFF-SCREEN BUFFERING)
% Decouples parallel math from serial frame dumping to prevent cluster stalls.
% =========================================================================

function generate_heatmap_videos_parallel(matFile, label, outDir, fps)
if nargin < 4 || isempty(fps), fps = 5; end
if ~exist(outDir, 'dir'), mkdir(outDir); end

stride = 1;

% Enforce pure software rendering across main process and all workers
warning('off', 'all');
try, opengl('save', 'software'); catch; end

pool = gcp('nocreate');
if isempty(pool)
    fprintf('Starting MATLAB parallel pool...\n');
    parpool(); 
end

currDir = fileparts(mfilename('fullpath'));
if isempty(currDir), currDir = pwd; end

pctRunOnAll(sprintf('cd(''%s'');', currDir));
pctRunOnAll(sprintf('addpath(''%s'', ''-begin'');', currDir));
pctRunOnAll clear functions;
pctRunOnAll rehash;
pctRunOnAll warning('off', 'all');

S = load(matFile);
if isfield(S, 'out')
    out = S.out;
elseif isfield(S, 'outRestart')
    out = S.outRestart;
else
    error('generate_heatmap_videos_parallel:noOut', 'No ''out'' or ''outRestart'' found in %s', matFile);
end

nSteps = out.stopStep;
if nSteps < 1
    error('generate_heatmap_videos_parallel:noSteps', '%s has stopStep=%d', matFile, nSteps);
end

fields = {'pressure', @plot_select_native2d_pressure; ...
          'stress',   @plot_select_native2d_stress; ...
          'velocity', @plot_select_native2d_velocity};

tempFrameDir = fullfile(outDir, ['temp_frames_' label]);
if ~exist(tempFrameDir, 'dir'), mkdir(tempFrameDir); end

sampledSteps = 1:stride:nSteps;
numSampled = numel(sampledSteps);

for f = 1:size(fields, 1)
    fieldName = fields{f,1};
    plotFn    = fields{f,2};
    fprintf('\n=== Rendering [%s] (1 per %d steps): %s ===\n', label, stride, fieldName);

    fieldTempDir = fullfile(tempFrameDir, fieldName);
    if ~exist(fieldTempDir, 'dir'), mkdir(fieldTempDir); end

    stateHist_sliced = out.stateHist(1:nSteps);
    fluidHist_sliced = out.fluidHist(1:nSteps);

    if isfield(out, 'PHist') && ~isempty(out.PHist)
        PHist_sliced = out.PHist(:,:,1:nSteps);
    else
        PHist_sliced = zeros(size(fluidHist_sliced{1}.P,1), size(fluidHist_sliced{1}.P,2), nSteps);
        for stIdx = 1:nSteps, PHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.P; end
    end

    if isfield(out, 'RPHist') && ~isempty(out.RPHist)
        RPHist_sliced = out.RPHist(:,:,1:nSteps);
    else
        RPHist_sliced = zeros(size(fluidHist_sliced{1}.meshF.Rp,1), size(fluidHist_sliced{1}.meshF.Rp,2), nSteps);
        for stIdx = 1:nSteps, RPHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.meshF.Rp; end
    end

    if isfield(out, 'ZPHist') && ~isempty(out.ZPHist)
        ZPHist_sliced = out.ZPHist(:,:,1:nSteps);
    else
        ZPHist_sliced = zeros(size(fluidHist_sliced{1}.meshF.Zp,1), size(fluidHist_sliced{1}.meshF.Zp,2), nSteps);
        for stIdx = 1:nSteps, ZPHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.meshF.Zp; end
    end

    if isfield(out, 'urCHist') && ~isempty(out.urCHist)
        urCHist_sliced = out.urCHist(:,:,1:nSteps);
    else
        urCHist_sliced = zeros(size(fluidHist_sliced{1}.urC,1), size(fluidHist_sliced{1}.urC,2), nSteps);
        for stIdx = 1:nSteps, urCHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.urC; end
    end

    if isfield(out, 'uzCHist') && ~isempty(out.uzCHist)
        uzCHist_sliced = out.uzCHist(:,:,1:nSteps);
    else
        uzCHist_sliced = zeros(size(fluidHist_sliced{1}.uzC,1), size(fluidHist_sliced{1}.uzC,2), nSteps);
        for stIdx = 1:nSteps, uzCHist_sliced(:,:,stIdx) = fluidHist_sliced{stIdx}.uzC; end
    end

    outBase = struct();
    outBase.meshE      = out.meshE;      
    outBase.meshL      = out.meshL;      
    outBase.par        = out.par;        
    outBase.t          = out.t;          
    outBase.stopStep   = nSteps;         
    outBase.fluidHist  = fluidHist_sliced;
    outBase.stateHist  = stateHist_sliced;
    outBase.PHist      = PHist_sliced;
    outBase.RPHist     = RPHist_sliced;
    outBase.ZPHist     = ZPHist_sliced;
    outBase.urCHist    = urCHist_sliced;
    outBase.uzCHist    = uzCHist_sliced;

    if isfield(out, 'cfg'),        outBase.cfg        = out.cfg; end
    if isfield(out, 'meshF'),      outBase.meshF      = out.meshF; end
    if isfield(out, 'interfaceE'), outBase.interfaceE = out.interfaceE; end
    if isfield(out, 'interfaceL'), outBase.interfaceL = out.interfaceL; end
    if isfield(out, 'z'),          outBase.z          = out.z; end
    if isfield(out, 'zGrid'),      outBase.zGrid      = out.zGrid; end
    if isfield(out, 'dz'),         outBase.dz         = out.dz; end
    if isfield(out, 'dtHist'),     outBase.dtHist     = out.dtHist; end

    % Single persistent figure for master process to prevent worker canvas locking
    fig = figure('Visible', 'off', 'Units', 'pixels', 'Position', [100 100 1200 800], ...
                 'GraphicsSmoothing', 'off', 'Renderer', 'painters');

    for frameIdx = 1:numSampled
        k = sampledSteps(frameIdx);
        clf(fig);
        ax = gca(fig);
        
        try
            outK = outBase;
            outK.stopStep = k;
            outK.state = stateHist_sliced{k};
            outK.fluid = fluidHist_sliced{k};

            if isfield(stateHist_sliced{k}, 'uE'), outK.state.uE = stateHist_sliced{k}.uE; end
            if isfield(stateHist_sliced{k}, 'uL'), outK.state.uL = stateHist_sliced{k}.uL; end

            if isfield(out, 'cfg') && isfield(out.cfg, 'fluid')
                outK.cfg.fluid = out.cfg.fluid;
            else
                outK.cfg.fluid.useFull2DFluid = true;
                outK.cfg.fluid.useHybridGap1DExterior2DFluid = false;
            end

            Pk   = PHist_sliced(:,:,k);
            Rk   = RPHist_sliced(:,:,k);
            Zk   = ZPHist_sliced(:,:,k);
            urk  = urCHist_sliced(:,:,k);
            uzk  = uzCHist_sliced(:,:,k);

            outK.native2D = struct();
            outK.native2D.P     = Pk;
            outK.native2D.R     = Rk;
            outK.native2D.Z     = Zk;
            outK.native2D.ur    = urk;
            outK.native2D.uz    = uzk;

            if strcmp(fieldName, 'stress')
                plotFn(outK, k, fig);
            else % pressure & velocity
                plotFn(outK, k, ax);
                set(ax, 'Units', 'normalized', 'Position', [0.12 0.12 0.75 0.80]);
                rMax = max(Rk(:)) * 1e6;
                if ~isfinite(rMax) || rMax <= 0, rMax = 6.0; end
                xlim(ax, [0 rMax]);
                ylim(ax, [-6 10]);
            end

            imgFile = fullfile(fieldTempDir, sprintf('frame_%06d.png', frameIdx));
            
            % Native off-screen vector printing (100% lockup proof on headless Linux)
            print(fig, imgFile, '-dpng', '-r100', '-opengl');

        catch MEk
            fprintf('Warning: Frame step %d (frame %d) skipped: %s\n', k, frameIdx, MEk.message);
        end

        if mod(frameIdx, 20) == 0 || frameIdx == numSampled
            fprintf('  Rendered %d / %d frames...\n', frameIdx, numSampled);
        end
    end

    if ishandle(fig), close(fig); end

    imgFiles = dir(fullfile(fieldTempDir, 'frame_*.png'));
    if isempty(imgFiles)
        fprintf('  [Warning] No frames rendered for %s; skipping video output.\n', fieldName);
        continue;
    end

    tempAviFile = fullfile(outDir, [label '_' fieldName '_temp.avi']);
    finalMp4File = fullfile(outDir, [label '_' fieldName '.mp4']);

    v = VideoWriter(tempAviFile, 'Uncompressed AVI');
    v.FrameRate = fps;
    open(v);

    firstImg = imread(fullfile(fieldTempDir, imgFiles(1).name));
    [targetH, targetW, ~] = size(firstImg);

    for idx = 1:numel(imgFiles)
        img = imread(fullfile(fieldTempDir, imgFiles(idx).name));
        if size(img, 1) ~= targetH || size(img, 2) ~= targetW
            img = imresize(img, [targetH, targetW]);
        end
        writeVideo(v, img);
    end
    close(v);

    % Conversion to MP4 using native MPEG-4 encoder
    cmd = sprintf(['ffmpeg -y -i "%s" ' ...
                   '-c:v mpeg4 -b:v 8M ' ...
                   '-pix_fmt yuv420p ' ...
                   '-vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" "%s"'], ...
        tempAviFile, finalMp4File);
        
    [status, cmdOut] = system(cmd);

    if status == 0 && exist(finalMp4File, 'file')
        delete(tempAviFile);
        fprintf('  Successfully compiled smooth MP4: %d frames -> %s\n', numel(imgFiles), finalMp4File);
    else
        actualAviFile = fullfile(outDir, [label '_' fieldName '.avi']);
        movefile(tempAviFile, actualAviFile);
        fprintf('  [Warning] FFmpeg conversion notice (%s). Saved AVI: %s\n', strtrim(cmdOut), actualAviFile);
    end
end

if exist(tempFrameDir, 'dir'), rmdir(tempFrameDir, 's'); end
fprintf('\nHeatmap video generation completed successfully for: %s\n', label);
end
