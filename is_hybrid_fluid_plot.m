function tf = is_hybrid_fluid_plot(out)
%IS_HYBRID_FLUID_PLOT
% Gate used by plot_select_native2d_stress.m, plot_select_native2d_pressure.m,
% plot_select_native2d_velocity.m, plot_final_native2d_pressure.m, and
% plot_final_native2d_velocity.m to decide whether there's a plottable
% native-2D fluid field to contour.
%
% BUG FIX: originally this only returned true for
% fluid.meshType == 'hybrid_gap1d_exterior2d'. But there are exactly two
% valid fluid mesh types in this codebase (see solve_fluid_2D_bodyfitted_MAC.m
% and solve_fluid_hybrid_gap1d_exterior2d_withnodes.m):
%   'hybrid_gap1d_exterior2d'  (the "1D" hybrid gap-1D/exterior-2D solve)
%   'bodyfitted_MAC'           (the "2D" full body-fitted MAC solve)
% Both produce a perfectly valid, plottable native-2D fluid field
% (meshF, ur2D, uz2D, pCell, etc.) -- so restricting this to only the
% hybrid type meant every one of the plotting functions above silently
% skipped drawing the fluid entirely whenever called on a full-2D
% ('bodyfitted_MAC') run, leaving a blank/white gap where the fluid
% contour should be. This is exactly the "there is no fluid inside!" gap
% seen in the 2D step-1 stress plot.
%
% Despite its name, this function's real job is "does out have a
% plottable native-2D fluid field," not literally "is this the hybrid
% solver" -- so it now accepts either valid mesh type.

    tf = false;
    if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
            ~isfield(out, 'fluidHist') || numel(out.fluidHist) < out.stopStep || ...
            isempty(out.fluidHist{out.stopStep})
        return;
    end
    fluidPlot = out.fluidHist{out.stopStep};
    validMeshTypes = {'hybrid_gap1d_exterior2d', 'bodyfitted_MAC'};
    tf = isfield(fluidPlot, 'meshType') && ...
        any(strcmpi(fluidPlot.meshType, validMeshTypes));
end
