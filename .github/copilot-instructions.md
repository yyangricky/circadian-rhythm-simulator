# Copilot / AI Agent Instructions

Purpose: Help AI coding agents quickly become productive in this repository by describing the project's structure, common patterns, run/debug workflows, and small gotchas discovered in the source.

- **Quick start**: Install the Python dependency and run the example scripts.
  - **Install**: `pip install circadian`
  - **Run examples**: `python actogram.py` and `python amplitude_plot.py`

- **Big picture / data flow**:
  - Input light schedules are generated (functions or `LightSchedule` classes) -> passed into model `equilibrate` to get initial conditions -> model integration called as `model(time, init, light_values)` or `model.integrate(...)` -> trajectories produce daily markers via `model.dlmos()` or `model.cbt(...)` -> plotting and metrics functions consume trajectories (`Actogram`, `plot_phase_amp`, `recovery_time`, etc.).

- **Key files to read first**:
  - `actogram.py` — quick demo that constructs `time`, `LightSchedule` instances, runs several models (`Forger99`, `Jewett99`, `Hannay19`, `Hannay19TP`), and plots an actogram.
  - `amplitude_plot.py` — the primary analysis harness: builds schedules, runs many models in `MODEL_NAMES`, computes metrics (re-entrainment, 90% amplitude), and calls `plot_phase_amp()`.
  - `recovery_time.py` — simple utility that finds first day where a DLMO series returns within a threshold of baseline.

- **Project-specific conventions & patterns**:
  - Units: time is always in **hours**. Arrays and constants like `DT` are fractional hours (e.g., `DT = 0.1`). Be careful converting minutes vs hours.
  - Baseline & recovery indexing: many scripts use `baseline_days = 30` and a recovery start at `start_day = 31` (or similar). Check constants like `HOURS`, `DAY_LUX`, `ALL_NIGHTER_*` when changing scenarios.
  - Tolerances: `TOL_MIN` is defined in minutes in `amplitude_plot.py` (e.g., `TOL_MIN = 15`) and converted to hours in functions that compare phase error (`tol_min / 60.0`).
  - Model imports: code expects an installed `circadian` package (`circadian.models`, `circadian.plots`, `circadian.lights`) but `amplitude_plot.py` falls back to `import models` if `circadian.models` isn't found. The local `circadian/` directory in the workspace is currently empty — verify whether you should use the PyPI package or local development copy.

- **Small gotchas & important notes**:
  - `recovery_time.py` docstring claims a default threshold of 10% but the function signature sets `threshold=0.30`. Confirm intended default before changing behavior.
  - `actogram.py` runs top-level code (no `if __name__ == "__main__"` guard). Edits that import `actogram.py` will execute side-effects — prefer adding a `main()` function and guard before refactors.
  - Plotting uses blocking `plt.show()`; for CI or automated runs prefer `plt.savefig()` or set headless backend `MATPLOTLIB_BACKEND=Agg`.

- **Where to modify behavior** (concrete examples):
  - Change time resolution: update `time = np.arange(0, 24*total_days, 0.10)` in `actogram.py` (change `0.10` to coarser/finer `DT`).
  - Swap or subset models: `actogram.py` and `amplitude_plot.py` instantiate models by class name; modify the list `MODEL_NAMES` in `amplitude_plot.py` or comment out model instances in `actogram.py`.
  - Adjust recovery criteria: `TOL_MIN`, `STREAK`, and `threshold` (in `recovery_time.py`) control re-entrainment detection.

- **Testing / debugging workflow**:
  - No unit tests are present. Run scripts directly to verify behavior.
  - Use prints/logging around these key calls: `model.equilibrate(...)`, `model(...)` integration, `model.dlmos()` or `model.cbt(...)` to inspect intermediate arrays.
  - If import fails for `circadian.models`, check whether a local `models.py` is intended (look for `import models` fallback in `amplitude_plot.py`).

- **Examples (copy-paste friendly)**:
  - Run the actogram demo: `python actogram.py`
  - Quick edit to reduce runtime (larger step): in `actogram.py` change `time = np.arange(0, 24*total_days, 0.10)` -> `time = np.arange(0, 24*total_days, 0.5)`
  - To avoid GUI in scripts (CI):
    - set env var before running: `SETX MATPLOTLIB_BACKEND Agg` (PowerShell: use `$Env:MATPLOTLIB_BACKEND = 'Agg'` for session)

- **If you are making structural changes**:
  - Preserve the current dataflow: light schedule -> equilibrate -> integrate -> extract markers/metrics -> plot. Refactor incremental functions rather than rewrite end-to-end scripts.
  - Add `if __name__ == "__main__": main()` guards before converting simple demo scripts into importable modules.

If anything here is unclear or you'd like the file to include more concrete code snippets (e.g., exact refactor for adding `main()` to `actogram.py`), tell me which section to expand and I'll iterate.
