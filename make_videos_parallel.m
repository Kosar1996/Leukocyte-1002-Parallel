%% =========================================================================
% AUDITED & ROBUST PARALLEL VIDEO LAUNCHER: make_videos_parallel.m
% =========================================================================

softlubeDir = fileparts(mfilename('fullpath'));
if isempty(softlubeDir)
    softlubeDir = pwd;
end

addpath(softlubeDir);

pctRunOnAll(sprintf('addpath(''%s'');', softlubeDir));
pctRunOnAll warning('off', 'MATLAB:graphics:noGraphicsAcceleration');

% Prioritize non-empty, valid dataset
inputname = fullfile(softlubeDir, 'case_0924_v2.mat');
if ~exist(inputname, 'file') || dir(inputname).bytes == 0
    inputname = fullfile(softlubeDir, 'case_0924_v2.mat');
end

if ~exist(inputname, 'file') || dir(inputname).bytes == 0
    matFiles = dir(fullfile(softlubeDir, '*.mat'));
    validMat = false;
    for m = 1:numel(matFiles)
        if matFiles(m).bytes > 1000
            inputname = fullfile(softlubeDir, matFiles(m).name);
            validMat = true;
            break;
        end
    end
    if ~validMat
        error('No valid non-empty .mat output file found in %s.', softlubeDir);
    end
end

fprintf('Target dataset selected: %s\n', inputname);

outputDir = fullfile(softlubeDir, 'videos');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

generate_heatmap_videos_parallel(inputname, '0924_', outputDir, 5);
