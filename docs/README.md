# Documentation

This folder contains the technical report for the HIF–iron homeostasis
simulation code.

## Build the PDF

Requires a LaTeX distribution (TeX Live, MiKTeX, MacTeX, …) with `pdflatex`.

```bash
cd docs
pdflatex Report.tex
pdflatex Report.tex
```

The output is `Report.pdf`.

## Contents

| Chapter | Topic |
|---------|--------|
| 1 | Overview and architecture |
| 2 | Installation and usage |
| 3 | Historical sampling protocol (lab notebook, cleaned) |
| 4 | Stiffness investigation notes |

The original long laboratory notes live in the sibling folder
`../Report_simulations/` (outside this cleaned package). This `docs/` version
is the maintained, publication-oriented documentation.
