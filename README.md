# Cellular HIF–Iron Homeostasis Network

Parameter-space exploration and ODE simulation for a cellular HIF–iron
homeostasis model.

## Authors

Affiliation at the time of this work: **TIMC-IMAG, CNRS**, Grenoble, France.

- **Eric Fanchon** and **Alexandre Rocca** — original simulation framework
  (sampling, ODE integration, constraint pipeline).
- **Irina Piazza** — software redesign, bug fixes, I/O/plotting utilities, and
  extension toward the HIF–iron model (this public version).

Code comments that start with `# IP <date>` mark Irina’s changes.

## Citation

DOI: *to be added*

## License

This project is released under the [MIT License](LICENSE).

## Contributions (Irina Piazza)

- Bug fixes in the simulation path (commented in place with `# IP …`).
- Introduced `config.py` for constants and settings previously scattered in the code.
- Clearer terminal display of parameters/values via pandas `DataFrame`.
- Introduced `read_file.py` for unambiguous bound / results input.
- Redesigned input/output around pandas DataFrames.
- Introduced `log.py` (run metadata and copy of inputs into each `Run_*` folder).
- Rewrote `filter.py` for a clearer post-filtering pipeline.
- Changed how the multiple iron replenishment experiments are passed and run.
- Unified the ODE right-hand side: one `diff_eqs` function with a
  `perturbation` flag (instead of two nearly identical steady-state / perturbed
  functions).
- Added `dump_solution_into_files()`, `compute_Fe_tot()`, and improved
  resampling / plotting workflow.
- Added plotting utilities: `plot_results.py`, `plot_params_correlation.py`,
  `multiplot_*_params_correlation.py`.

## Installation

### 1. Prerequisites

- **Python 3.8 or newer** (3.10+ recommended)
- `pip` and `venv` (included with standard Python installs)

Check your Python version:

```bash
python --version
```

On some systems the command is `python3` instead of `python`.

### 2. Create a virtual environment

From the project root (`hif-iron-homeostasis`):

**Linux / macOS**

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
pip install -r requirements.txt
```

**Windows (PowerShell)**

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
pip install -r requirements.txt
```

**Windows (Command Prompt)**

```bat
python -m venv .venv
.venv\Scripts\activate.bat
python -m pip install --upgrade pip
pip install -r requirements.txt
```

To leave the environment later: `deactivate`.

### 3. Optional: Conda

```bash
conda create -n hif-iron python=3.11
conda activate hif-iron
pip install -r requirements.txt
```

### 4. Verify the install

```bash
python tests/run_all_tests.py
```

All tests should pass (`OK`).

## Project layout

| File / folder | Role |
|---------------|------|
| `simulate.py` | Main entry point (sampling, simulation, filtering, plots) |
| `model.py` | ODE system and biological constraints |
| `config.py` | Model constants, bounds, solver settings |
| `filter.py` | Post-filtering of accepted solutions |
| `read_file.py` | Bound / results file I/O |
| `log.py` | Run metadata (`info.txt`) |
| `pca.py`, `sobol_lib.py` | Sampling helpers |
| `plot_*.py`, `multiplot_*.py` | Offline plotting utilities |
| `b-vars-hif.txt`, `b-par-hif*.txt` | Variable / parameter bound boxes |
| `tests/` | Contract + regression tests |
| `docs/` | Extended documentation (LaTeX report) |

## Run a simulation

With the virtual environment activated:

```bash
python simulate.py b-vars-hif.txt b-par-hif.txt results.txt Box-Sampling 1000
```

Arguments:

1. variable bounds file  
2. parameter bounds file (or previous results, depending on mode)  
3. output results file name  
4. mode (`Box-Sampling`, `Re-Sampling`, `PCA-guided-sampling`, `Shell-Sampling`, …)  
5. number of points  

Each run writes a timestamped `Run_*` folder with solutions, plots, and `info.txt`.

Filter is called automatically after sampling when solutions are found.

## Tests

```bash
python tests/run_all_tests.py
```

Regression goldens live in `tests/golden/`. Refresh them only after an
intentional numerical change:

```bash
python tests/generate_golden.py
```

## Documentation

A more detailed technical report (architecture, workflows, historical sampling
protocol, stiffness notes) is in [`docs/`](docs/).

Build the PDF (requires a LaTeX distribution such as TeX Live or MiKTeX):

```bash
cd docs
pdflatex Report.tex
pdflatex Report.tex
```
