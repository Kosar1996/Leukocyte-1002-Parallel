# Changes 1008

- `monolithic_two_solids_residual_jacobian_scaled.m`: the active-force activation fix (provided 10/7; see the serial
  repository) plus Set 4: the active-force tangent is assembled with one sparse() call from triplets instead of
  16 indexed sparse updates per element (same element values and order; bit-identical, unit test test_active_1008).
- `run_0930_par.m`, `run_0930_prof.m`: optional POOL_TYPE=threads (thread-based pool, Set 3; default unchanged:
  process pool) and optional test switches SOFTLUBE_ACTIVE_FORCE, SOFTLUBE_F0, SOFTLUBE_ZBOT (unset = values in the file).
- Previous versions: `pre_change_backup_1008/`.
