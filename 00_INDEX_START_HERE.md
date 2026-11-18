# CpG Island Detector: File Index (start here)

Python 3.8+ , standard library only. Tested on Python 3.12.3 (77 tests, all OK).
Outputs are computationally predicted candidate regions, not confirmed CpG islands.

## Labelled file list
| Label | File | What it is |
|---|---|---|
| DOC-0 | `00_INDEX_START_HERE.md` | this index |
| DOC-1 | `README.md` | installation, usage, options, exit codes |
| DOC-2 | `report/01_project_report_science_and_design.md` | overview, scientific foundation, maths, algorithm and alternatives, validation, limitations, references |
| DOC-3 | `report/02_code_walkthrough.md` | function-by-function explanation of the code, with hand-traced examples |
| CODE-1 | `src/cpg_detector.py` | the program (importable module and command-line tool) |
| CODE-2 | `scripts/generate_sample_fasta.py` | reproducible generator for the synthetic sample |
| TEST-1 | `tests/test_cpg_detector.py` | 77 unittest tests with hand-calculated expected values |
| DATA-1 | `data/sample.fasta` | 6 SYNTHETIC records with known island positions |
| CONFIG-1 | `requirements.txt` | states there are no third-party dependencies |
| RESULT-1 | `results/sample_run_output.txt` | real detector output on DATA-1 |
| RESULT-2 | `results/test_run_log.txt` | real log of the test run |
| RESULT-3 | `results/sample_regions.csv` | CSV export of predicted regions |
| RESULT-4 | `results/plots/001...006_*.svg` | SVG profiles, one per sequence (open in a browser) |

## Map to the requested deliverables
1. Overview and objectives: DOC-2 section 1
2. Scientific foundation and maths: DOC-2 sections 2-3
3. Algorithm and alternatives: DOC-2 section 4
4. Directory structure: DOC-1
5. Source code: CODE-1, CODE-2 (explained in DOC-3)
6. Sample FASTA: DATA-1
7. Test suite: TEST-1
8. Install and run: DOC-1
9. Example session and output: RESULT-1, RESULT-2
10. Limitations and future work: DOC-2 sections 7-8

## Quick start
    python -m unittest discover -s tests -v
    python src/cpg_detector.py -f data/sample.fasta -o results/regions.csv --plot-dir results/plots
