# Code Walkthrough: `src/cpg_detector.py`

This document explains the program function by function, so that you can
understand, trace and defend it even without running it. Real captured outputs
are in `results/sample_run_output.txt` (detector on `data/sample.fasta`) and
`results/test_run_log.txt` (77 tests, all OK, Python 3.12.3).

## 1. Big picture: data flow

```
CLI args ──► load_records ──► validate_sequence / parse_fasta_lines
                                  │ (SequenceRecord: id, description, clean sequence)
                                  ▼
                           analyse_record
        ┌─────────────────────────┼──────────────────────────────┐
        ▼                         ▼                              ▼
  compute_metrics      PrefixCounts + generate_window_bounds   evaluate_criteria
  (whole sequence)     (score every window in O(1))            (pass / reasons)
                                  │ qualifying windows
                                  ▼
                           merge_intervals ─► re-measure ─► CandidateRegion (1-based)
                                  ▼
          format_report · write_csv / write_json · build_svg_profile
```

## 2. Module header and constants
* `from __future__ import annotations` lets type hints work on Python 3.8.
* Defaults `DEFAULT_MIN_LENGTH_BP=200`, `DEFAULT_MIN_GC_PERCENT=50.0`,
  `DEFAULT_MIN_OBS_EXP=0.6` are the historical Gardiner-Garden & Frommer values.
  `DEFAULT_WINDOW_SIZE_BP=200` and `DEFAULT_STEP_SIZE_BP=10` are practical scan
  settings, not biological constants.
* `UNAMBIGUOUS_BASES`, `IUPAC_AMBIGUOUS_BASES`, `ALL_ACCEPTED_SYMBOLS` define
  which characters are legal. All other numbers (plot sizes, limits) are named
  constants too, so there are no magic numbers in function bodies.

## 3. Exceptions
`CpGDetectorError` is the base for deliberate user-facing errors.
`InvalidSequenceError`, `FastaFormatError` and `ParameterError` (also a
`ValueError`) subclass it. `main()` catches the base class and prints
`error: <message>` with exit code 1, so users never see a traceback for bad input.

## 4. Data structures (frozen dataclasses)
* **`SequenceRecord(identifier, description, sequence)`**: one validated,
  upper-case sequence.
* **`DetectionParameters`**: criteria and scan settings. `__post_init__`
  validates every value (e.g. step ≥ 1, GC between 0 and 100) and raises
  `ParameterError`, so an invalid configuration cannot exist.
* **`SequenceMetrics`**: raw counts (`length`, `a/c/g/t_count`,
  `ambiguous_count`, `cpg_count`) plus derived properties:
  * `unambiguous_length` = A+C+G+T (the "L" in the formulas).
  * `gc_percent` = 100·(C+G)/L, `None` if L = 0.
  * `obs_exp_cpg` = CpG·L/(C·G), `None` if C = 0 or G = 0.
  * `ambiguous_fraction` = ambiguous/length.
  Returning `None` (not 0, not an exception) makes "undefined" explicit; later
  code treats `None` as failing.
* **`CandidateRegion(record_id, start, end, metrics)`**: 1-based inclusive
  coordinates; `length = end - start + 1`.
* **`WindowPoint`**: per-window GC/Obs/Exp used only for plotting.
* **`RecordResult`**: bundles the record, whole-sequence metrics, the reasons the
  whole sequence fails (if any), the regions and the optional profile.

## 5. Coordinates: `to_one_based_inclusive(start0, end0_exclusive)`
Returns `(start0 + 1, end0_exclusive)`. Example: Python slice `[400, 800)` →
reported `401–800`. Reverse: `sequence[start-1:end]`. This is the only place the
conversion happens, so conventions cannot drift.

## 6. Input handling

### `validate_sequence(raw, allow_ambiguous=True, label="sequence")`
1. Remove all whitespace (`"".join(raw.split())`) and upper-case. Lower-case is
   a legal representation of the same nucleotide, so this is not a silent
   "repair" of an invalid base.
2. `set(cleaned) - ALL_ACCEPTED_SYMBOLS` finds illegal characters quickly (O(n),
   no Python loop on valid input).
3. If any exist, build a message listing each as `'X' at position 5` (first 10
   only, with total count) plus a hint for `U` (RNA). Raises
   `InvalidSequenceError`. Nothing is ever converted to A/C/G/T.
4. If ambiguity codes are disallowed, raise a separate error naming the first
   position.
Positions refer to the sequence *after* whitespace removal.

### `parse_fasta_lines(lines, source, allow_ambiguous)`
A small state machine: `header` holds the current record's header, `chunks` its
sequence lines.
* Blank lines are skipped; `\r\n` is handled by `strip()`.
* A line starting with `>` closes the previous record via the inner
  `finish_record()` and opens a new one. Empty header → `FastaFormatError`.
* Sequence text before any header → error citing the line number.
* `finish_record()` splits the header into identifier (first word) and
  description, rejects a record with no sequence, validates the joined sequence
  (error message names the record) and appends it.
* No records at all (empty file) → error.

### `read_fasta_file(path, allow_ambiguous)`
Opens with `utf-8-sig` (tolerates a BOM). A binary/non-text file raises
`UnicodeDecodeError`, converted to a clear `FastaFormatError`. A missing file
raises `FileNotFoundError`, which `main()` reports as `file not found: <path>`.

## 7. Sequence analysis

### `compute_metrics(sequence)`
Validates (so invalid characters can never be hidden inside the "ambiguous"
count), then uses `str.count` for A, C, G, T and `"CG"`.
`ambiguous_count = length − (A+C+G+T)`. `str.count("CG")` is exactly the number
of adjacent C→G pairs because "CG" cannot overlap itself: `"CGCG"` → 2,
`"CGCGCG"` → 3. `"GC"` is not counted, and `"CNG"` has no CpG (N sits between).

**Worked example `AACGCGTT`:** A=2, C=2, G=2, T=2, CpG=2 →
GC = 4/8 = 50 %, Obs/Exp = 2·8/(2·2) = 4.0. This was run:
`GC 50.00 %, CpG 2, Obs/Exp 4.000`, and it fails the default criteria only
because GC is not *strictly* above 50 % (and the sequence is shorter than 200).

### `find_cpg_positions(sequence)`
Repeated `str.find("CG", i+1)`; returns 1-based positions of each C. For
`AACGCGTT` → `[3, 5]`. Used for the CpG-density strip in the plot.

### `evaluate_criteria(metrics, params, check_length=True)`
Returns a *list of human-readable failure reasons*; an empty list means pass.
Checks: length ≥ minimum (skippable for windows, which may be shorter than the
minimum island length); GC strictly greater; Obs/Exp strictly greater; undefined
values fail; ambiguous fraction ≤ limit. Returning reasons (not just a boolean)
lets the report explain *why* a sequence fails, and lets tests assert that
exactly one criterion failed.

## 8. Region detection

### `PrefixCounts(sequence)`
Builds five cumulative-count arrays (`array('i')`, 4 bytes per entry) for A, C,
G, T and CpG in one pass. `cpg[k]` = number of CpGs whose C has 0-based index < k;
a C counts as a CpG start if the next base is G.
`metrics(start, end)` for the slice `[start, end)`:
* base counts = `prefix[end] − prefix[start]`;
* CpGs fully inside have their C at indices `start … end−2`, so
  `cpg[end−1] − cpg[start]` (0 if the window is shorter than 2).
A CpG straddling a window boundary is therefore counted in neither window (tested
on `ACGT`). One test checks that this O(1) calculation equals direct
`compute_metrics` on the slice for many windows of a pseudo-random sequence with
N, i.e. two independent implementations must agree.

### `generate_window_bounds(sequence_length, params)`
Generator of 0-based end-exclusive `(start, end)` pairs.
* Shorter than `min_length` → no windows (nothing can qualify).
* Not longer than the window → one window covering everything.
* Else starts 0, step, 2·step, … ≤ n−W. If the last grid start is not n−W, one
  extra window flush with the end is yielded, so the 3′ tail is never skipped.
Example n=450, W=200, S=100 → `(0,200) (100,300) (200,400) (250,450)`.

### `merge_intervals(intervals)`
Sort, then merge when `next_start <= current_end` (overlapping **or touching**).
`[(300,500),(0,100),(50,150),(150,200),(400,450)]` → `[(0,200),(300,500)]`.

### `analyse_record(record, params, include_profile=False)`
1. `compute_metrics` + `evaluate_criteria` on the whole sequence (for the
   report's "Whole sequence meets all criteria" line).
2. For each window, get O(1) metrics; if it passes (length check off), remember
   it. Optionally record a `WindowPoint` for plotting.
3. Merge qualifying windows.
4. **Re-measure each merged region as a whole and drop it unless it passes every
   criterion including length.** The union of passing windows does not
   automatically pass.
5. Convert to 1-based inclusive and build `CandidateRegion`s.

Hand-traced example from the tests: 400 bp of `AT`, 400 bp of `CGCGCGAT`, 400 bp
of `AT`; W=200, S=100. Window [300,500) holds 100 island bases → GC 37.5 % →
fails. Windows [400,600), [500,700), [600,800) lie inside the island → pass.
[700,900) fails again. The three passing windows overlap and merge to [400,800)
→ reported **401–800**, length 400, 150 CpG, GC 75 %, Obs/Exp = 150·400/(150·150) = 2.667.

`detect_regions(record_id, sequence, params)` is a convenience wrapper.

## 9. Reporting and export
* `render_text_map`: one-line map; each character covers `length/width` bases
  and is `#` if it overlaps any region.
* `format_report`: builds the text report (no timestamps → deterministic).
  `NA` is printed for undefined values.
* `write_csv` / `write_json` / `export_results`: CSV has one row per region with
  counts and the parameters used (self-describing); a `.json` path gives JSON
  with parameters, per-record summary and a disclaimer. Parent folders are
  created; newlines are fixed to `\n`.

## 10. Visualisation (`build_svg_profile`)
Dependency-free SVG with: GC% curve with the 50 % threshold dashed; Obs/Exp curve
(axis capped at 1.5) with the 0.6 threshold dashed; shaded predicted regions;
a region track; a CpG-density strip (300 bins, opacity proportional to count);
position axis. `None` values break the line instead of being drawn as zero. Long
sequences are downsampled to ≤ 1200 points. Record IDs are XML-escaped.
`safe_file_stem` turns identifiers into safe file names.
(Verified as well-formed XML by a test; not rendered to an image here, so check
in a browser.)

## 11. Command line (`build_parser`, `load_records`, `main`)
* Exactly one input source is required: `-s/--sequence`, `-f/--fasta` or
  `-i/--interactive`.
* Criteria, scan settings, `--ambiguous {allow,reject}`, `-o`, `--plot-dir`.
* `main()` builds `DetectionParameters` first; invalid values call
  `parser.error` → usage message, **exit 2**. Runtime/input problems print
  `error: …` to stderr, **exit 1**. Success is **exit 0**.
* An empty `-s ""` is allowed: the report shows `GC NA %` and 0 regions.

## 12. How to verify locally (if you cannot run it here)
```bash
python --version                      # 3.8 or newer
python -m unittest discover -s tests -v   # expect: Ran 77 tests ... OK
python src/cpg_detector.py -f data/sample.fasta
```
Compare against `results/sample_run_output.txt`. On the sample you should see 5
candidate regions in total; both negative controls and the short record give 0.
If your Python is 3.8–3.11, note that I only ran it on 3.12.3.
