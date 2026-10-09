# Changes 1009

- Ramped active force (provided 10/8): `softlube_run_case_global_coupled.m` scales the active force density with a
  raised-cosine ramp, f0_eff = f0 * 0.5 (1 - cos(pi tau)), tau = min(1, (t - t_start)/t_ramp), t_ramp = 3e-6 s;
  `monolithic_two_solids_residual_jacobian_scaled.m` and `solve_monolithic_two_solids_fsolve_timestep.m` no longer
  overwrite the ramped value (parL.fz_active_translocation is set from par only if parL does not have it).
- Ramp start (1009): t_start = the start time of the run (t = 0 for a fresh start, nargin == 1; the time of the restart
  step for a restart, nargin >= 2), instead of the hard-coded 3.25e-5 s (step 13). The log prints the ramp start and
  the ramp factor of every step ([1009 RAMP] lines, print only).
- Removed (10/8 request): the `SOFTLUBE_STOP_AFTER_STEP` environment-variable stop in
  `softlube_run_case_global_coupled.m` (section [1002]) and its mention in README.md. The simulation length is set
  only by nSteps in the run file.
- Previous versions: `pre_change_backup_1009/`.

Parallel repository: the same changes on top of Set 5 (the jacobian file keeps the Set 4 active-force assembly; only the one line above changes).
