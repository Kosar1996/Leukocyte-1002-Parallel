function report = check_interface_traction_mismatch_report(out, plotstep, opts)
%CHECK_INTERFACE_TRACTION_MISMATCH_REPORT

% traction on the two sides of the interface (leukocyte|fluid and
% fluid|endothelium), rotate both into the local normal/tangential frame
% of the interface curve r = delta(z), and report the mismatch
%   en(s) = tn_solid(s) - tn_fluid(s)
%   et(s) = tt_solid(s) - tt_fluid(s)
% as a function of interface arc length s, along with max/median/90th
% percentile of |en| and |et|.
%
% This reuses the exact rotation/interpolation machinery already in
% check_interface_traction_continuity.m (same stress recovery, same
% viscoelastic solid stress via recover_nodal_stress_axisym_viscoelastic.m,
% same local-slope normal-vector construction), but:
%   1. evaluates on a DENSE z-grid instead of 5 spot checks, so percentile
%      statistics are meaningful;
%   2. converts z to physical arc length s along each side's own deformed
%      interface curve;
%   3. returns numeric arrays + summary stats instead of only printing;
%   4. optionally plots en(s), et(s) for both sides.
%
% Usage:
%   report = check_interface_traction_mismatch_report(out, plotstep)
%   report = check_interface_traction_mismatch_report(out, plotstep, opts)
%
% opts (all optional):
%   opts.nQuery   number of interface sample points (default 61)
%   opts.epsFrac  fraction of local gap width used to offset the sample
%                 point off the wall on each side (default 0.1, same
%                 default as check_interface_traction_continuity.m)
%   opts.trimFrac fraction of the z-range trimmed off each end before
%                 sampling, to stay inside the fluid mesh's interpolation
%                 domain (default 0.05)
%   opts.makePlot true/false (default true)
%   opts.label    string used in plot titles / printed header, e.g.
%                 '1D hybrid, step 1' (default '')
%
% report fields (each a struct with sub-structs .leuko and .endo):
%   report.leuko.s / .en / .et         arc length [m] and mismatch [Pa] arrays
%   report.leuko.stats.maxAbsEn / .medianAbsEn / .p90AbsEn
%   report.leuko.stats.maxAbsEt / .medianAbsEt / .p90AbsEt
%   report.endo.*                      same fields for the endothelium side
%   report.t, report.plotstep, report.label

if nargin < 3 || isempty(opts), opts = struct(); end
if ~isfield(opts,'nQuery')   || isempty(opts.nQuery),   opts.nQuery   = 61;   end
if ~isfield(opts,'epsFrac')  || isempty(opts.epsFrac),  opts.epsFrac  = 0.1;  end
if ~isfield(opts,'trimFrac') || isempty(opts.trimFrac), opts.trimFrac = 0.05; end
if ~isfield(opts,'makePlot') || isempty(opts.makePlot), opts.makePlot = true; end
if ~isfield(opts,'label'),    opts.label = '';                                end
% Same rule-of-thumb "slowly varying" bound used in check_gap1DWindow_validity.m
% (|dr/dz| << 1). Exposed here so callers can flag which sample points sit in
% a "steep slope" region of the interface (where the lubrication approximation
% is least valid) vs a "low slope" region, instead of guessing from z alone.
if ~isfield(opts,'slopeThreshold') || isempty(opts.slopeThreshold)
    opts.slopeThreshold = 0.15;
end

fluid = out.fluidHist{plotstep};

st = out.stateHist{plotstep};
if ~isfield(st, 'uEPrev') || ~isfield(st, 'uLPrev')
    error(['check_interface_traction_mismatch_report: out.stateHist{%d} is missing ', ...
        'uEPrev/uLPrev -- cannot compute the viscous stress contribution.'], plotstep);
end
if isfield(out, 'dtHist') && numel(out.dtHist) >= plotstep && isfinite(out.dtHist(plotstep))
    dtStep = out.dtHist(plotstep);
else
    dtStep = out.par.dt;
    warning('out.dtHist not available for step %d; falling back to par.dt.', plotstep);
end

parL = out.par;
if isfield(out.par, 'GL')  && isfinite(out.par.GL),  parL.Ge = out.par.GL;  end
if isfield(out.par, 'KL')  && isfinite(out.par.KL),  parL.Ke = out.par.KL;  end
if isfield(out.par, 'etaL') && isfinite(out.par.etaL), parL.etaE = out.par.etaL; end

mismatchOpts = struct('nQuery', opts.nQuery, 'epsFrac', opts.epsFrac, 'trimFrac', opts.trimFrac);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, out.par, ...
    out.meshL, st.uL, st.uLPrev, parL, ...
    dtStep, fluid, mismatchOpts);

sL = leuko.s; enL = leuko.en; etL = leuko.et; slopeLv = leuko.slope;
sE = endo.s;  enE = endo.en;  etE = endo.et;

report = struct();
report.plotstep = plotstep;
report.label = opts.label;
if isfield(out,'t') && numel(out.t) >= plotstep
    report.t = out.t(plotstep);
else
    report.t = NaN;
end

report.leuko = leuko;
report.endo  = endo;

fprintf('\n=== Interface traction mismatch report: %s step %d (t=%.4g s) ===\n', ...
    opts.label, plotstep, report.t);
print_side('Leukocyte | fluid', report.leuko);
print_side('fluid | Endothelium', report.endo);

if opts.makePlot
    figure('Name', sprintf('Interface traction mismatch %s step %d', opts.label, plotstep));
    subplot(2,1,1);
    plot(sL*1e6, enL, 'o-', 'DisplayName', 'e_n leukocyte side'); hold on;
    plot(sE*1e6, enE, 's-', 'DisplayName', 'e_n endothelium side');
    yline(0, 'k:');
    shade_slope_regions(gca, sL*1e6, slopeLv, opts.slopeThreshold);
    xlabel('interface arc length s [\mum]'); ylabel('e_n = t_n^{solid} - t_n^{fluid}  [Pa]');
    legend('Location','best'); grid on;
    title(sprintf('Normal traction mismatch, %s step %d (t=%.4g s)', opts.label, plotstep, report.t));

    subplot(2,1,2);
    plot(sL*1e6, etL, 'o-', 'DisplayName', 'e_t leukocyte side'); hold on;
    plot(sE*1e6, etE, 's-', 'DisplayName', 'e_t endothelium side');
    yline(0, 'k:');
    shade_slope_regions(gca, sL*1e6, slopeLv, opts.slopeThreshold);
    xlabel('interface arc length s [\mum]'); ylabel('e_t = t_t^{solid} - t_t^{fluid}  [Pa]');
    legend('Location','best'); grid on;
    title(sprintf('Tangential traction mismatch (shaded = |dr/dz| > %.2f, "steep slope")', opts.slopeThreshold));
end

end

function shade_slope_regions(ax, sVals, slope, thresh)
% Shades contiguous "steep slope" bands (|dr/dz| > thresh, where the 1D
% lubrication approximation is least valid) behind the already-plotted
% data, and labels each band, so region membership doesn't have to be
% guessed from z alone. Must be called AFTER the real data is plotted, so
% ylim(ax) reflects the actual data range instead of MATLAB's default
% empty-axes [0,1] range.
yl = ylim(ax);   % capture the data-driven range before adding patches
steep = abs(slope(:)) > thresh;
d = diff([0; steep(:); 0]);
startIdx = find(d==1);
endIdx = find(d==-1) - 1;

for k = 1:numel(startIdx)
    sLo = sVals(startIdx(k));
    sHi = sVals(endIdx(k));
    p = patch(ax, [sLo sHi sHi sLo], [yl(1) yl(1) yl(2) yl(2)], ...
        [1 0.85 0.85], 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
        'HandleVisibility', 'off');
    uistack(p, 'bottom');   % keep shading behind the data lines/markers
    text(ax, 0.5*(sLo+sHi), yl(2), 'steep slope', ...
        'Color', [0.6 0 0], 'FontSize', 8, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top');
end
ylim(ax, yl);   % patches span exactly [yl(1) yl(2)], so this is a no-op
                % safety net against MATLAB re-autoscaling to include them
end

function print_side(name, side)
fprintf('  %-20s |en|: max=%9.4f  median=%9.4f  p90=%9.4f  Pa   (peak at r=%.3f, z=%.3f um)\n', ...
    name, side.stats.maxAbsEn, side.stats.medianAbsEn, side.stats.p90AbsEn, ...
    side.stats.locMaxEn(1)*1e6, side.stats.locMaxEn(2)*1e6);
fprintf('  %-20s |et|: max=%9.4f  median=%9.4f  p90=%9.4f  Pa   (peak at r=%.3f, z=%.3f um)\n', ...
    '', side.stats.maxAbsEt, side.stats.medianAbsEt, side.stats.p90AbsEt, ...
    side.stats.locMaxEt(1)*1e6, side.stats.locMaxEt(2)*1e6);
end
