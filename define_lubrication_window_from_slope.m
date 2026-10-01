function [gapZ, info] = define_lubrication_window_from_slope( ...
    meshL, uL, interfaceL, meshE, uE, interfaceE, opts)
%DEFINE_LUBRICATION_WINDOW_FROM_SLOPE
% Computes the 1D lubrication window [zLo, zHi] (i.e. gap1DWindow /
% cfg.fluid.hybridGapZ) from the INSTANTANEOUS deformed leukocyte and
% endothelium interface geometry, instead of a fixed constant chosen once
% and never revisited.
%
% BACKGROUND: the 1D Reynolds/lubrication equation used inside the gap1D
% window assumes the gap height h(z) = deltaE(z) - deltaL(z) varies slowly
% with z, i.e. |dr/dz| << 1 on BOTH bounding interfaces. plot_hybrid_vs_2D_
% final_comparison.m and check_gap1DWindow_validity.m established that the
% fixed default gap1DWindow = [-0.2e-6, 4.2e-6] extends into z-ranges where
% the leukocyte interface slope exceeds that assumption, and that
% restricting the window to where the assumption actually holds reduced
% the leukocyte-interface traction mismatch by 40-58%. This function
% automates picking that window from the current deformed shape, so it can
% be recomputed at any step instead of hand-tuned once.
%
% METHOD:
%   1. Sample r(z) on both deformed interfaces over their z-overlap, on a
%      shared fine query grid (so the two sides -- which may have
%      different native node spacing -- are compared apples-to-apples).
%   2. Compute dr/dz on that shared grid for both sides.
%   3. A point is "safe" (lubrication-valid) only if BOTH sides satisfy
%      |dr/dz| < opts.slopeThreshold -- violating either bounding surface's
%      slow-variation assumption breaks the approximation.
%   4. Find all contiguous safe z-ranges ("candidates").
%   5. Choose the candidate containing opts.referenceZ if one exists,
%      otherwise the widest candidate overall.
%
% Usage:
%   gapZ = define_lubrication_window_from_slope(meshL, uL, interfaceL, meshE, uE, interfaceE)
%   [gapZ, info] = define_lubrication_window_from_slope(..., opts)
%
%   % Recompute from a solved run at any step k:
%   gapZ = define_lubrication_window_from_slope( ...
%       out.meshL, out.stateHist{k}.uL, out.interfaceL, ...
%       out.meshE, out.stateHist{k}.uE, out.interfaceE);
%
% opts (all optional):
%   opts.slopeThreshold  |dr/dz| bound for "slowly varying" (default 0.15,
%                        the same rule-of-thumb threshold used in
%                        check_gap1DWindow_validity.m and
%                        check_interface_traction_mismatch_report.m).
%   opts.referenceZ      z-value [m] the chosen window should contain if
%                        possible -- typically the current minimum-gap
%                        location, i.e. the physically active part of the
%                        interaction. Default 0.
%   opts.nQuery          number of points on the shared query grid used to
%                        evaluate slope (default 401).
%   opts.minWidth        minimum acceptable window width [m]; if the
%                        chosen candidate is narrower, a warning is issued
%                        (the window is still returned as-is; the caller
%                        decides how to handle a too-narrow result).
%                        Default 0 (no minimum enforced).
%
% Outputs:
%   gapZ  1x2 [zLo, zHi] in meters. Errors if no safe range exists at all
%         (i.e. the lubrication approximation is not valid ANYWHERE in the
%         current interface-overlap region).
%   info  struct with diagnostic detail:
%           .candidates     Nx2 list of all safe windows found, sorted by
%                            zLo
%           .chosenIdx      row of .candidates that was returned as gapZ
%           .zQuery         shared query grid used
%           .drdzL, .drdzE  slope on the query grid, both sides
%           .safe           logical mask, true where both sides are safe

if nargin < 7 || isempty(opts), opts = struct(); end
if ~isfield(opts,'slopeThreshold') || isempty(opts.slopeThreshold)
    opts.slopeThreshold = 0.15;
end
if ~isfield(opts,'referenceZ') || isempty(opts.referenceZ)
    opts.referenceZ = 0;
end
if ~isfield(opts,'nQuery') || isempty(opts.nQuery)
    opts.nQuery = 401;
end
if ~isfield(opts,'minWidth') || isempty(opts.minWidth)
    opts.minWidth = 0;
end

[zL, rL] = deformed_interface_rz(meshL, uL, interfaceL);
[zE, rE] = deformed_interface_rz(meshE, uE, interfaceE);

zLo = max(min(zL), min(zE));
zHi = min(max(zL), max(zE));
if ~(zHi > zLo)
    error(['define_lubrication_window_from_slope: leukocyte and endothelium ', ...
        'interfaces do not overlap in z (leukocyte z=[%.4g, %.4g], ', ...
        'endothelium z=[%.4g, %.4g]).'], min(zL), max(zL), min(zE), max(zE));
end

zQuery = linspace(zLo, zHi, opts.nQuery);
rLq = interp1(zL, rL, zQuery, 'linear');
rEq = interp1(zE, rE, zQuery, 'linear');
drdzL = gradient(rLq, zQuery);
drdzE = gradient(rEq, zQuery);

safe = abs(drdzL) < opts.slopeThreshold & abs(drdzE) < opts.slopeThreshold;

candidates = safe_ranges(zQuery, safe);

info = struct();
info.candidates = candidates;
info.zQuery = zQuery(:);
info.drdzL = drdzL(:);
info.drdzE = drdzE(:);
info.safe = safe(:);

if isempty(candidates)
    error(['define_lubrication_window_from_slope: the lubrication ', ...
        'approximation (|dr/dz| < %.3g on both interfaces) is not valid ', ...
        'anywhere in the current z-overlap [%.4g, %.4g] m. No window to return.'], ...
        opts.slopeThreshold, zLo, zHi);
end

containsRef = candidates(:,1) <= opts.referenceZ & opts.referenceZ <= candidates(:,2);
if any(containsRef)
    chosenIdx = find(containsRef, 1, 'first');
else
    widths = candidates(:,2) - candidates(:,1);
    [~, chosenIdx] = max(widths);
end

gapZ = candidates(chosenIdx, :);
info.chosenIdx = chosenIdx;

if (gapZ(2) - gapZ(1)) < opts.minWidth
    warning(['define_lubrication_window_from_slope: chosen window width ', ...
        '%.4g m is narrower than opts.minWidth = %.4g m.'], ...
        gapZ(2) - gapZ(1), opts.minWidth);
end

end

function [z, r] = deformed_interface_rz(mesh, u, interfaceIds)
ids = interfaceIds(:);
r = mesh.nodes(ids,1) + u(2*ids-1);
z = mesh.nodes(ids,2) + u(2*ids);
[z, ord] = sort(z);
r = r(ord);
end

function ranges = safe_ranges(z, safe)
d = diff([0; safe(:); 0]);
startIdx = find(d==1);
endIdx = find(d==-1) - 1;
ranges = [z(startIdx), z(endIdx)];
end
