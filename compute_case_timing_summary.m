function row = compute_case_timing_summary(matFile, caseLabel, jobElapsedSeconds, nProcessors)
%COMPUTE_CASE_TIMING_SUMMARY Build one summary-table row for item 4
% (simulation time, average time/timestep, number of processors).
%
% Usage:
%   row = compute_case_timing_summary('/path/out_case_i.mat', '(i) Base', 7397, 1)
%
% matFile           : path to a saved case output (.mat with 'out' or 'outRestart')
% caseLabel          : label for the table row, e.g. '(i) Base'
% jobElapsedSeconds  : total wall-clock seconds for the Slurm job (from sacct
%                       Elapsed), or [] / NaN if not yet known
% nProcessors        : number of CPUs used for this job
%
% Average time/timestep excludes the final accepted step (per instruction)
% and any step whose wall time is a statistical outlier vs. the rest
% (median * 3 threshold), flagged as "stuck". If the run predates the
% stepWallTimeHist instrumentation, falls back to
% jobElapsedSeconds / nSteps as an approximation and flags it as such.

S = load(matFile);
if isfield(S, 'out')
    out = S.out;
elseif isfield(S, 'outRestart')
    out = S.outRestart;
else
    error('compute_case_timing_summary:noOut', 'No ''out'' or ''outRestart'' variable found in %s', matFile);
end

nSteps = out.stopStep;

row = struct();
row.case = caseLabel;
row.nProcessors = nProcessors;
row.nStepsCompleted = nSteps;

if nargin < 3 || isempty(jobElapsedSeconds) || ~isfinite(jobElapsedSeconds)
    row.simTimeSeconds = NaN;
    row.simTimeNote = 'pending -- fill in from sacct Elapsed once job completes';
else
    row.simTimeSeconds = jobElapsedSeconds;
    row.simTimeNote = '';
end

hasStepTiming = isfield(out, 'stepWallTimeHist') && numel(out.stepWallTimeHist) >= nSteps && nSteps >= 1;

if hasStepTiming && nSteps >= 2
    t = out.stepWallTimeHist(1:nSteps);
    t = t(:);
    usable = t(1:end-1);   % exclude final step per instruction

    med = median(usable(~isnan(usable)));
    isStuck = usable > 3 * med;
    nStuckExcluded = sum(isStuck);
    usableClean = usable(~isStuck);

    row.avgTimePerStep = mean(usableClean);
    row.avgTimePerStepNote = sprintf('from stepWallTimeHist, excluded final step + %d stuck-step outlier(s)', nStuckExcluded);
elseif isfinite(row.simTimeSeconds) && nSteps >= 1
    row.avgTimePerStep = row.simTimeSeconds / nSteps;
    row.avgTimePerStepNote = 'APPROXIMATION: run predates per-step timing instrumentation; total elapsed / nSteps (includes final step, cannot exclude stuck steps)';
else
    row.avgTimePerStep = NaN;
    row.avgTimePerStepNote = 'pending -- need either stepWallTimeHist (new solver) or sacct Elapsed time';
end

fprintf('%-25s procs=%-3d steps=%-4d simTime=%s avgTime/step=%s (%s)\n', ...
    row.case, row.nProcessors, row.nStepsCompleted, ...
    ternary_str(isfinite(row.simTimeSeconds), sprintf('%.1fs', row.simTimeSeconds), 'pending'), ...
    ternary_str(isfinite(row.avgTimePerStep), sprintf('%.2fs', row.avgTimePerStep), 'pending'), ...
    row.avgTimePerStepNote);
end

function s = ternary_str(cond, a, b)
if cond
    s = a;
else
    s = b;
end
end
