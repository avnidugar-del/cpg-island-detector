# CpG Island Detector

A Python command-line tool that predicts **candidate CpG-rich regions** in DNA
sequences using the historical Gardiner-Garden & Frommer (1987) criteria
(length, GC content, observed/expected CpG ratio) with a sliding-window scan
and window merging.

> **Scientific honesty:** outputs are *computationally predicted candidate
> regions* based on sequence composition only. They are **not** experimentally
> confirmed CpG islands and say nothing by themselves about promoter function,
> methylation status, gene expression or disease. See `report/01_project_report_science_and_design.md`.

## Requirements
Python 3.8+ and nothing else (standard library only).
Developed and tested with Python 3.12.3.

## Layout
```
cpg-island-detector/
    README.md
    requirements.txt
    data/sample.fasta                  # SYNTHETIC sequences with known structure
    scripts/generate_sample_fasta.py   # reproducible generator for the sample
    src/cpg_detector.py                # the whole tool (importable + CLI)
    tests/test_cpg_detector.py         # unittest suite (77 tests)
    results/                           # outputs (includes captured sample run + test log)
    report/01_project_report_science_and_design.md  # science, algorithm, validation, limits
    report/02_code_walkthrough.md         # function-by-function explanation of the code
```

## Usage
```bash
# FASTA file (multiple records OK), CSV export and SVG plots
python src/cpg_detector.py -f data/sample.fasta -o results/regions.csv --plot-dir results/plots

# Sequence typed on the command line / at a prompt
python src/cpg_detector.py -s CGCGCGATCGCGCGAT --min-length 10 --window-size 10 --step-size 2
python src/cpg_detector.py -i

# Stricter, Takai-Jones-style criteria
python src/cpg_detector.py -f data/sample.fasta --min-length 500 --min-gc 55 --min-obs-exp 0.65 --window-size 500

# JSON export (chosen by the .json extension)
python src/cpg_detector.py -f data/sample.fasta -o results/regions.json
```
Run `python src/cpg_detector.py --help` for all options.

| Option | Default | Meaning |
|---|---|---|
| `--min-length` | 200 | minimum region length (bp, `>=`) |
| `--min-gc` | 50 | GC % must be strictly greater |
| `--min-obs-exp` | 0.6 | Obs/Exp CpG must be strictly greater |
| `--window-size` / `--step-size` | 200 / 10 | scan settings (practical, not biological, values) |
| `--max-ambiguous-fraction` | 0.05 | windows/regions with more N etc. are ineligible |
| `--ambiguous` | allow | `reject` makes any IUPAC ambiguity code an error |

Exit codes: 0 success, 1 input/runtime error (message on stderr), 2 bad
command-line usage or invalid parameter.

## Tests
```bash
python -m unittest discover -s tests -v
```

## Coordinates
All output is **1-based, inclusive** (`401-800` = bases 401..800, length 400).
Internally Python slices are 0-based, end-exclusive: `sequence[start-1:end]`.

## Interpreting outputs
* **Whole sequence** line: metrics over the entire input, and whether the whole
  thing meets the criteria (with reasons if not).
* **Candidate regions**: merged qualifying windows, re-measured as a whole and
  re-checked against *all* criteria. Boundaries are accurate only to roughly
  the step size and can overshoot a true boundary by up to about a window.
* **Obs/Exp = NA** means undefined (no C or no G); it never counts as passing.
* **SVG plot**: GC% and Obs/Exp curves with dashed threshold lines, shaded
  predicted regions, and a CpG-density strip.
