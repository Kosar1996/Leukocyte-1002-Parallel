# Leukocyte-1005-Parallel

Base: Leukocyte_Main_Files-1001_v1 (stable active force) + the same changes as Leukocyte-1005-Serial + parfor.
Originals of all changed files: `pre_change_backup_1005/`. Results are bit-identical to the serial code for any number of workers.

| File | Change |
|---|---|
| `assemble_finite_def_axisym.m`, `assemble_finite_def_internal_force_only.m`, `assemble_axisym_kelvin_voigt_viscous.m`, `assemble_axisym_kelvin_voigt_viscous_force_only.m` | `parfor` over blocks of consecutive elements (`[1005 OPT]`, one block per worker; `SOFTLUBE_PARFOR_BLOCKS` overrides), serial accumulation in the original element order; `assemble_finite_def_axisym` with one output skips the element tangent and `K` |
| `solve_finite_def_solid.m` | `[1005 OPT]` line-search trial residuals with one output (same file as the serial repo) |
| `softlube_run_case_global_coupled.m` | Restart fix `[1002 FIX]` on her v1 solver (same file as the serial repo) |
| `apply_interface_traction.m`, `apply_interface_traction_sensitivity.m` | Faster traction (same files as the serial repo; parfor off by default, `TRACTION_PARFOR=1` to switch on) |
| `run_0930_par.m` | `run_0930.m` + pool block (`NUM_PROCS` workers, own pool folder per Slurm job, client `maxNumCompThreads(2)`) + run options `RESTART_FILE`, `RESTART_STEP`, `FRESH_START=1` |
| `run_0930_prof.m` | `run_0930_par.m` with the MATLAB profiler (Master Table 4) |

Run: `NUM_PROCS=8` in the environment, then `run('run_0930_par.m')`. Input `.mat` files are not in the repo.
