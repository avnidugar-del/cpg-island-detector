# CpG Island Detector

Two dependency-free Python tools that predict **candidate CpG islands** in DNA from sequence composition alone, using the classic Gardiner-Garden & Frommer (1987) criteria: length ≥ 200 bp, GC > 50 %, observed/expected CpG > 0.6. Both scan with sliding windows, merge the passing windows, and re-check each merged region against every criterion.

> Outputs are computationally predicted candidates, **not** experimentally confirmed CpG islands. They say nothing about methylation, promoter activity, expression or disease. All sample data are simulated or constructed, and neither tool has been benchmarked on real genomes.

## What's inside

| | Implementation A | Implementation B |
|---|---|---|
| Folder | `implementation-a/` | `implementation-b/` |
| Highlights | strict validation, SVG profile plots | presets, boundary trimming, BED export |
| Tests | 77 | 29 |
| Outputs | text, CSV, JSON, SVG | console, TSV, BED, JSON |

Python 3.8+, standard library only (developed and run on 3.12.3).

## Run it

```bash
# A: tests, then detect with CSV and SVG plots
(cd implementation-a && python -m unittest discover -s tests)
python implementation-a/src/cpg_detector.py -f implementation-a/data/sample.fasta \
    -o regions.csv --plot-dir plots

# B: tests, then detect with TSV, BED and JSON
(cd implementation-b && python -m unittest discover -s tests)
python implementation-b/cpg_detector.py implementation-b/data/demo.fasta \
    --tsv out.tsv --bed out.bed --json out.json
```

## Example

On a 600-bp test record (100 bp flank, 400 bp island, 100 bp flank) the whole record has GC 46.67 % and fails the 50 % criterion, yet scanning in windows finds the island: `101-498` (B) and `51-550` (A, untrimmed, so it overhangs the designed `101-500`).

## Learn more

* Science, maths and algorithm: `implementation-a/report/` and `implementation-b/docs/report.md`
* Each program run on the other's data: `docs/CROSS_VALIDATION.md`
* Limitations: thresholds are conventions, fixed windows blur boundaries, forward strand only, no repeat masking, and the references in the reports were written from memory and need checking before citing.
