%% CHECK_WEAK_FORCE_BALANCE
% The radial/axial debug run PROVED the traction-correction converges to a
% genuine displacement fixed point (|dRadial|,|dAxial| -> exact 0), while
% the POINTWISE %mismatch (EnE etc.) stays stuck at ~78%. That means the
% correction is enforcing one notion of "solid matches fluid" (finite-
% element force balance under the applied traction load) while the
% mismatch check uses a DIFFERENT notion (pointwise agreement between two
% independently-recovered stress fields, sampled at query points).
%
% This checks whether the two sides are actually close in the NET/
% INTEGRATED (bulk) sense -- i.e. does the total force from the solid's
% recovered traction roughly balance the total force from the fluid's
% recovered traction along the interface, even though the pointwise
% comparison looks bad? If bulk force balance is much better than the
% pointwise 78%, that points at a stress-RECOVERY/measurement artifact
% (e.g. noisy nodal extrapolation near the traction boundary) rather than
% genuinely unbalanced physics -- exactly the "measuring the wrong thing"

%
% Reuses compute_interface_traction_mismatch.m's own per-point solid/fluid
% traction values (endo.en = tnSolid-tnFluid, endo.tnFluid = tnFluid, so
% tnSolid = endo.en + endo.tnFluid) -- same recovery/interpolation
% pipeline already validated elsewhere, just integrated instead of
% compared pointwise. Same for leukocyte (leuko) and tangential (et/tt).

clc; close all;
cd(fileparts(mfilename('fullpath')));

cfg = build_cfg_full2D_pressure2('nSteps', 1, 'useHybridGap1DExterior2DFluid', true);
testRelax = 0.3;
testMaxIter = 15;
cfg.fluid.bodyFittedTractionCorrectionRelax = testRelax;
cfg.parOverrides.bodyFittedTractionCorrectionRelax = testRelax;
cfg.fluid.maxBodyFittedTractionCorrections = testMaxIter;
cfg.parOverrides.maxBodyFittedTractionCorrections = testMaxIter;

fprintf('Running relax=%.2f, maxIterations=%d to reach the converged fixed point ...\n', testRelax, testMaxIter);
out = softlube_run_case_global_coupled(cfg);

st = out.stateHist{1};
fl = out.fluidHist{1};
par = out.par;

parLmismatch = par;
if isfield(par, 'GL') && isfinite(par.GL), parLmismatch.Ge = par.GL; end
if isfield(par, 'KL') && isfinite(par.KL), parLmismatch.Ke = par.KL; end
if isfield(par, 'etaL') && isfinite(par.etaL), parLmismatch.etaE = par.etaL; end

mismatchOpts = struct('nQuery', 61, 'epsFrac', 0.1, 'trimFrac', 0.05);
[leuko, endo] = compute_interface_traction_mismatch( ...
    out.meshE, st.uE, st.uEPrev, par, out.meshL, st.uL, st.uLPrev, parLmismatch, par.dt, fl, mismatchOpts);

report_side('Endothelium', endo);
report_side('Leukocyte', leuko);

fprintf('\nSMOKE_TEST_STATUS: OK\n');

function report_side(name, side)
tnSolid = side.en + side.tnFluid;
ttSolid = side.et + side.ttFluid;

maskN = isfinite(tnSolid) & isfinite(side.tnFluid);
maskT = isfinite(ttSolid) & isfinite(side.ttFluid);
nDroppedN = numel(tnSolid) - nnz(maskN);
nDroppedT = numel(ttSolid) - nnz(maskT);

netNormalSolid = trapz(side.s(maskN), tnSolid(maskN));
netNormalFluid = trapz(side.s(maskN), side.tnFluid(maskN));
netTangentSolid = trapz(side.s(maskT), ttSolid(maskT));
netTangentFluid = trapz(side.s(maskT), side.ttFluid(maskT));

refN = max(abs(trapz(side.s(maskN), abs(side.tnFluid(maskN)))), 1e-9);
refT = max(abs(trapz(side.s(maskT), abs(side.ttFluid(maskT)))), 1e-9);

pctNetNormal = 100 * abs(netNormalSolid - netNormalFluid) / refN;
pctNetTangent = 100 * abs(netTangentSolid - netTangentFluid) / refT;

fprintf('\n=== %s ===\n', name);
fprintf('  Pointwise (p90):  normal=%.2f%%  tangential=%.2f%%  (dropped %d/%d normal, %d/%d tangential NaN query points)\n', ...
    side.stats.pctP90En, side.stats.pctP90Et, nDroppedN, numel(tnSolid), nDroppedT, numel(ttSolid));
fprintf('  Net/integrated force along interface:\n');
fprintf('    normal:      solid=%.6e  fluid=%.6e  |diff|/|fluid|=%.2f%%\n', ...
    netNormalSolid, netNormalFluid, pctNetNormal);
fprintf('    tangential:  solid=%.6e  fluid=%.6e  |diff|/|fluid|=%.2f%%\n', ...
    netTangentSolid, netTangentFluid, pctNetTangent);
end
