# Changes 1008

- `monolithic_two_solids_residual_jacobian_scaled.m`: the active-force activation fix (provided 10/7; see the serial
  repository) plus Set 4: the active-force tangent is assembled with one sparse() call from triplets instead of
  16 indexed sparse updates per element (same element values and order; bit-identical, unit test test_active_1008).
- `run_0930_par.m`, `run_0930_prof.m`: optional POOL_TYPE=threads (thread-based pool, Set 3; default unchanged:
  process pool) and optional test switches SOFTLUBE_ACTIVE_FORCE, SOFTLUBE_F0, SOFTLUBE_ZBOT (unset = values in the file).
- Previous versions: `pre_change_backup_1008/`.

## Set 5 (same day, on the profiled bottleneck)

- `solve_finite_def_solid.m`, `assemble_finite_def_axisym.m`: the line search of the solid Newton solver (75% of
  step 14 with the active force) evaluates its trials (alpha = 1, 1/2, 1/4, ...) one after the other, each in its own
  short parallel call. Set 5 assembles the residuals of up to 10 consecutive trials in ONE parallel call and then
  checks them in the original order with the original acceptance test, so the accepted step is the same and the
  results are bit-identical (unit test test_set5_1008; runs N = 1-16 bit-identical to the reference). Trials after
  the accepted one are computed but not used. SOFTLUBE_LS_BATCH=1 restores the original loop; without a parallel
  pool the original loop is used.
- Previous versions of these two files: `pre_change_backup_1008/`.
