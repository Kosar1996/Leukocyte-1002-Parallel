function diag = check_pressure_stress_consistency(out, plotstep)
%CHECK_PRESSURE_STRESS_CONSISTENCY  Way-1 diagnostic.
% Quantifies how well the recovered fluid stress trace,
% -(1/3)*trace(sigma), matches the solver's own pressure field p, and
% separately reports the solver's own mass-conservation residual
% (divInf), so a mismatch between p and -(1/3)*trace(sigma) can be
% attributed to the post-processing stress recovery step rather than to
% an inaccurate flow solve.
%
% Background: for an incompressible Newtonian fluid,
%   sigma = -p*I + 2*mu*D  =>  trace(sigma) = -3p + 2*mu*div(u)
% so -(1/3)*trace(sigma) = p only where div(u) = 0 exactly. The MAC
% solver enforces mass conservation on its native staggered grid; the
% stress recovery in recover_fluid_nodes_pressure_stress_Q4.m instead
% differentiates an already-averaged cell-center velocity field
% (urC/uzC = 0.5*(face values), see solve_stokes_bodyfitted_MAC.m lines
% ~295-296) on a derived dual mesh. That is a different discrete
% divergence operator, so a small residual between p and
% -(1/3)*trace(sigma) is expected even when the flow solve itself is
% essentially exact.
%
%   diag = check_pressure_stress_consistency(out, plotstep)
%
% Fields of diag:
%   diffMax        max_e |pCell(e) - (-(1/3)*trace(sigmaCell(e,:)))| [Pa]
%   maxP           max(|pCell|)                                      [Pa]
%   diffRelToP     diffMax / maxP
%   viscousScale   representative magnitude of the normal viscous
%                  stress components (|sigma_ii - (-p)| = |2*mu*rate|),
%                  used to judge diffMax against something other than
%                  the overall pressure scale
%   diffRelToVisc  diffMax / viscousScale
%   divInf         solver's own max |div(u)| residual                [1/s]
%   divInfAsPa     (2*mu/3) * divInf, converted to the same units as
%                  diffMax for a direct side-by-side comparison
%   mu             viscosity used [Pa s]

fluid = out.fluidHist{plotstep};
mu = out.par.mu;

% Reuse the stored recovery if present, otherwise recompute it exactly
% the way softlube_run_case_global_coupled.m / read_and_plot.m do.
meshF = add_fluid_nodes(fluid.meshF);
[~, sigmaCell, center, pRaw] = recover_fluid_nodes_pressure_stress_Q4( ...
    meshF, fluid.ur2D, fluid.uz2D, mu, fluid.pCell);

traceThird = -(sigmaCell(:,1) + sigmaCell(:,2) + sigmaCell(:,3)) / 3;
diffAll = abs(pRaw(:) - traceThird(:));
valid = isfinite(diffAll);

diag = struct();
diag.diffMax = max(diffAll(valid));
diag.maxP = max(abs(pRaw(valid)));
diag.diffRelToP = diag.diffMax / max(diag.maxP, eps);

% Representative viscous normal-stress scale: |sigma_ii - (-p)|,
% i.e. |2*mu*(that normal strain rate)|, median over finite, nonzero
% cells so a handful of near-zero cells don't dominate the ratio.
viscNormal = abs(sigmaCell(:,1:3) + pRaw(:));
viscNormalFlat = viscNormal(:);
%viscNormalFinite = viscNormalFlat(isfinite(viscNormalFlat) & viscNormalFlat > 0);
viscNormalFlat = viscNormal(:);
viscNormalFinite = viscNormalFlat(isfinite(viscNormalFlat) & viscNormalFlat > 0);
if isempty(viscNormalFinite)
    diag.viscousScale = NaN;
else
    diag.viscousScale = median(viscNormalFinite);
end
diag.diffRelToVisc = diag.diffMax / max(diag.viscousScale, eps);
[~, idxMax] = max(diffAll(valid));
validIdx = find(valid);
idxMax = validIdx(idxMax);

diag.rrViscousAtMax = sigmaCell(idxMax,1) + pRaw(idxMax);   % 2*mu*durdr at that point
diag.ttViscousAtMax = sigmaCell(idxMax,2) + pRaw(idxMax);   % 2*mu*(ur/r) at that point
diag.diffRelToRRAtMax = diag.diffMax / max(abs(diag.rrViscousAtMax), eps);
diag.diffRelToTTAtMax = diag.diffMax / max(abs(diag.ttViscousAtMax), eps);
% Where is the global max actually happening?
diag.maxLocation_r_um = center(idxMax,1) * 1e6;
diag.maxLocation_z_um = center(idxMax,2) * 1e6;
fprintf('  location of max mismatch (r, z)     = (%.3f, %.3f) um\n', ...
    diag.maxLocation_r_um, diag.maxLocation_z_um);

% Restrict to the 1D lubrication gap window only -- the leukocyte/
% endothelium interaction region, not the coarser exterior blocks.
gapZ = [-0.2e-6, 4.2e-6];
inGap = center(:,2) >= gapZ(1) & center(:,2) <= gapZ(2) & valid;
if any(inGap)
    gapIdx = find(inGap);
    [diag.diffMaxInGap, kGap] = max(diffAll(inGap));
    idxGapMax = gapIdx(kGap);
    diag.rrViscousInGapAtMax = sigmaCell(idxGapMax,1) + pRaw(idxGapMax);
    diag.ttViscousInGapAtMax = sigmaCell(idxGapMax,2) + pRaw(idxGapMax);
    fprintf('  max mismatch WITHIN gap z=[%.2f,%.2f] um = %.6e Pa (at r=%.3f, z=%.3f um)\n', ...
        gapZ(1)*1e6, gapZ(2)*1e6, diag.diffMaxInGap, ...
        center(idxGapMax,1)*1e6, center(idxGapMax,2)*1e6);
    fprintf('  local viscous sigma_rr / sigma_tt in gap = %.6e / %.6e Pa\n', ...
        diag.rrViscousInGapAtMax, diag.ttViscousInGapAtMax);
    fprintf('  gap mismatch as %% of local rr / tt viscous = %.4f %% / %.4f %%\n', ...
        100*diag.diffMaxInGap/max(abs(diag.rrViscousInGapAtMax),eps), ...
        100*diag.diffMaxInGap/max(abs(diag.ttViscousInGapAtMax),eps));
end

% The worst single cell is landing at a corner where the gap/exterior
% stitching seam meets the outer domain wall -- likely an
% interpolation artifact, not representative of the interior. Look at
% the distribution, and at the max with boundary/seam-adjacent cells
% excluded, for a fairer comparison.
diag.diffPercentiles = prctile(diffAll(valid), [50 75 90 95 99]);
fprintf('  mismatch percentiles [50 75 90 95 99]  = ');
fprintf('%.4e  ', diag.diffPercentiles); fprintf('Pa\n');

rMarginUm = 0.3;   % exclude cells within this many microns of r = rOuter
zMarginUm = 0.3;   % exclude cells within this many microns of the gap-window edges
rOuterUm = max(center(:,1)) * 1e6;
interior = valid & ...
    (rOuterUm - center(:,1)*1e6 > rMarginUm) & ...
    (center(:,2)*1e6 - gapZ(1)*1e6 > zMarginUm) & ...
    (gapZ(2)*1e6 - center(:,2)*1e6 > zMarginUm);
if any(interior)
    diag.diffMaxInterior = max(diffAll(interior));
    fprintf('  max mismatch excluding cells within %.2f um of domain/seam edges = %.6e Pa\n', ...
        rMarginUm, diag.diffMaxInterior);
end

fprintf('  viscous sigma_rr AT the max-mismatch point = %.6e Pa\n', diag.rrViscousAtMax);
fprintf('  viscous sigma_tt AT the max-mismatch point = %.6e Pa\n', diag.ttViscousAtMax);
fprintf('  mismatch as %% of local rr viscous stress   = %.4f %%\n', 100*diag.diffRelToRRAtMax);
fprintf('  mismatch as %% of local tt viscous stress   = %.4f %%\n', 100*diag.diffRelToTTAtMax);

if isfield(fluid, 'divInf') && isfinite(fluid.divInf)
    diag.divInf = fluid.divInf;
else
    diag.divInf = NaN;
end
diag.divInfAsPa = (2*mu/3) * diag.divInf;
diag.mu = mu;

fprintf('\n--- Pressure / stress-trace consistency check (step %d) ---\n', plotstep);
fprintf('  max|p - (-1/3 trace(sigma))|        = %.6e Pa\n', diag.diffMax);
fprintf('  max|p|                              = %.6e Pa\n', diag.maxP);
fprintf('  ratio to max|p|                     = %.4f %%\n', 100*diag.diffRelToP);
fprintf('  representative viscous normal stress = %.6e Pa\n', diag.viscousScale);
fprintf('  ratio to viscous normal stress       = %.4f %%\n', 100*diag.diffRelToVisc);
fprintf('  solver max|div(u)| (divInf)         = %.6e 1/s\n', diag.divInf);
fprintf('  divInf converted to Pa (2*mu/3*div) = %.6e Pa\n', diag.divInfAsPa);
if isfinite(diag.divInf) && diag.divInfAsPa < 1e-3 * diag.diffMax
    fprintf(['  --> solver residual is >>1000x smaller than the observed ', ...
        'mismatch: the flow solve itself is essentially divergence-free, ', ...
        'so the mismatch is coming from the post-processing stress ', ...
        'recovery step, not the solve.\n']);
end
fprintf('------------------------------------------------------------\n\n');

end
