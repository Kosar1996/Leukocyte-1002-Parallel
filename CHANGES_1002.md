# Leukocyte-1002-Parallel

Base: Leukocyte_Main_Files-1001 (active force in the Jacobian solve) + the same fixes as Leukocyte-1002-Serial + parfor.
Originals: `pre_fix_backup_1002/` (solver, traction), `pre_parallelization_backup_1002/` (4 assembly files).

| File | Change |
|---|---|
| `assemble_finite_def_axisym.m`, `assemble_finite_def_internal_force_only.m`, `assemble_axisym_kelvin_voigt_viscous.m`, `assemble_axisym_kelvin_voigt_viscous_force_only.m` | `parfor` over elements, serial accumulation in the original order (bit-identical results for any number of workers) |
| `softlube_run_case_global_coupled.m` | Restart fix `[1002 FIX]` (same file as the serial repo) |
| `apply_interface_traction.m`, `apply_interface_traction_sensitivity.m` | Faster traction (same files as the serial repo; parfor off by default, `TRACTION_PARFOR=1` to switch on) |
| `run_0930_par.m` | `run_0930.m` + pool block (`NUM_PROCS` workers, own pool folder per Slurm job, client `maxNumCompThreads(2)`) + run options `RESTART_FILE`, `RESTART_STEP`, `FRESH_START=1` |

Run: `NUM_PROCS=8` in the environment, then `run('run_0930_par.m')`. Input `.mat` files are not in the repo.

## Update 10/2: restart state fix (`softlube_run_case_global_coupled.m`, marked `[1002 FIX]`)

A restart now continues exactly like the uninterrupted run, also after a step that needed a dt retry:
- next dt = min(par.dt, dtGrowFactor * dt of the restart step) after a retried/reduced step (was always par.dt);
- dtPrev of the fsolve warm start = dt of the restart step (was par.dt);
- a file written by the fixed code keeps the predictor setting (`useMonoPredictor`) of the run that wrote it; older files (e.g. the step-13 input) keep the original behaviour (predictor on).

Tests (check_fix_1002, OVERALL PASS): restart 13 -> 14 -> 15 bit-identical to 13 -> 15 (parallel N = 8 and serial); with the divergence retry (test code only, not in this repository) 0 -> 1 -> 2 bit-identical to 0 -> 2.
