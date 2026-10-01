function out = add_native2d_field(out)
%ADD_NATIVE2D_FIELD  Attach the out.native2D field that
% plot_select_native2d_stress.m (via draw_select_fluid_boundaries.m)
% needs, but which is only ever built automatically inside
% run_input_full2D_pressure2.m (out.native2D = collect_final_native2d(out)).
%
% When you call softlube_run_case_global_coupled(cfg) directly (e.g. via
% build_cfg_full2D_pressure2.m), out never gets this field, and
% plot_select_native2d_stress errors with "Unrecognized field name
% 'native2D'". This just replays that same construction here.
%
% NOTE: native2D always reflects the FINAL accepted step (out.stopStep) --
% it only supplies boundary/axis-limit context for the plot, not the
% actual per-step stress data. So call this once after the run finishes;
% you can then call plot_select_native2d_stress(out, k) for ANY step k,
% not just the last one.

if ~isfield(out, 'stopStep') || out.stopStep < 1 || ...
        ~isfield(out, 'PHist') || isempty(out.PHist)
    warning('add_native2d_field: out.PHist missing/empty; native2D not added.');
    return;
end

k = out.stopStep;
out.native2D.P     = out.PHist(:,:,k);
out.native2D.R     = out.RPHist(:,:,k);
out.native2D.Z     = out.ZPHist(:,:,k);
out.native2D.ur    = out.urCHist(:,:,k);
out.native2D.uz    = out.uzCHist(:,:,k);
out.native2D.speed = out.speedCHist(:,:,k);

end
