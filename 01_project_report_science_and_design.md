# CpG Island Detector: A Python-Based Computational Tool for Identifying CpG-Rich Regions in DNA Sequences

## 1. Overview and objectives

**Problem.** In vertebrate genomes, the dinucleotide 5'-CG-3' is much rarer than
the base composition predicts, except in short stretches called CpG islands.
These stretches are GC-rich, CpG-rich and usually unmethylated, and many sit at
gene promoters. Finding them from sequence alone is a classic computational
biology exercise.

**Objectives.**
1. Explain the biology and statistics behind CpG islands.
2. Implement a transparent, reproducible detector based on documented criteria.
3. Make every design decision explicit and defensible (thresholds, windows,
   ambiguous bases, coordinates).
4. Validate the code with independently calculated test values and with a
   synthetic sequence of known structure.
5. State plainly what the tool cannot tell you.

**What the tool is.** An educational/exploratory predictor of *candidate*
CpG-rich regions from sequence composition. **What it is not.** A tool that
finds experimentally confirmed CpG islands, or that says anything about
promoter activity, methylation, expression or disease.

---

## 2. Scientific foundation

### 2.1 CpG dinucleotides and the "p"
"CpG" means a cytosine immediately followed by a guanine along **one strand**,
read 5'→3'. The "p" stands for the **phosphate** (a phosphodiester bond)
linking the 3' carbon of the cytosine's sugar to the 5' carbon of the guanine's
sugar. It distinguishes this *consecutive-on-one-strand* arrangement from a
**C·G base pair**, where cytosine on one strand hydrogen-bonds to guanine on the
opposite strand. "CpG" ≠ "GpC": the latter is G followed by C.

CpG is its own reverse complement: the opposite strand reads 5'-CG-3' at the
same position. So the CpG count is the same on both strands, and the
methylation state of a CpG can be copied symmetrically after replication.

### 2.2 What CpG islands are
CpG islands are regions, typically a few hundred to a few thousand bp, that
differ from bulk DNA in three ways: higher GC content, CpG frequency close to
what base composition predicts (little depletion), and—in most somatic cells—
low DNA methylation. Bird and colleagues first described them in the 1980s
(initially as "HTF islands"). Gardiner-Garden & Frommer (1987) gave the
computational definition used here.

### 2.3 Biological relevance
* **Promoters.** A large share of vertebrate gene promoters (commonly cited as
  roughly 60–70 % in humans) lie in or near CpG islands, especially promoters of
  housekeeping and many developmental genes.
* **Methylation and gene regulation.** DNA methyltransferases add a methyl group
  to cytosine in CpG, forming 5-methylcytosine (5mC). Promoter-CpG-island
  methylation is generally associated with stable transcriptional silencing
  (e.g. X-chromosome inactivation, imprinting). Unmethylated islands are
  compatible with, but do not guarantee, transcriptional activity.
* **Epigenetics and disease.** In many cancers, specific promoter CpG islands
  (including those of some tumour-suppressor genes) become hypermethylated and
  silenced, while bulk genome methylation is reduced. Aberrant methylation is
  also studied in imprinting disorders and ageing.
* **Islands without promoters.** Many CpG islands are "orphan" or intragenic.
  Not every CpG-rich sequence is a regulatory element.

### 2.4 Why CpG is depleted outside islands
In bulk vertebrate DNA, most CpGs are methylated. 5mC can deaminate
spontaneously to **thymine**, producing a T·G mismatch. (Deamination of
unmethylated cytosine gives uracil, which cells recognise and remove
efficiently; a T·G mismatch is harder to recognise as damage.) Over evolutionary
time, such C→T changes accumulate in methylated CpGs, so CpG ends up at roughly a
fifth of its expected frequency in bulk DNA. CpG islands escape this because they
are protected from methylation, so their CpGs are retained.

### 2.5 Why GC content, Obs/Exp and length all matter
* **GC content** separates GC-rich islands from AT-rich bulk DNA.
* **Obs/Exp CpG** asks whether CpG is depleted *given* the amount of C and G.
  This is the signature that distinguishes an island from a merely GC-rich
  region: a GC-rich isochore can have high GC and still strong CpG depletion.
* **Length** matters statistically and biologically. In very short sequences,
  counts of C, G and CpG are tiny and the ratio is noisy; and very short CpG-rich
  stretches are not what is meant by an "island".

### 2.6 Why many CG dinucleotides are not enough
* **Short sequences** can show many CpGs by chance or by repetition (a 100-bp
  "CGCG…" repeat has Obs/Exp = 2.0 but is not 200 bp).
* **Obs/Exp can be high while GC is low** when C and G are rare but adjacent.
  The test `"CGATAT" × 34` has Obs/Exp = 6.0 and GC = 33 %.
* **CpG-rich repeats** (e.g. Alu elements, which are ≈300 bp, GC-rich and
  CpG-dense but typically methylated) can satisfy the 1987 criteria without being
  functional islands. This is why later work (Takai & Jones 2002) proposed
  stricter thresholds (≥ 500 bp, GC ≥ 55 %, Obs/Exp ≥ 0.65).

### 2.7 Limits of composition-only identification
Sequence composition cannot show methylation state, chromatin context,
promoter use, expression or disease relevance. Thresholds are conventions tuned
to human-like genomes; other genomes (GC-rich genomes, invertebrates with little
or no CpG methylation, plants with CHG/CHH methylation) behave differently.
Experimental island maps (e.g. affinity purification of unmethylated CpG
clusters, bisulfite sequencing) and annotations are needed for biological
conclusions.

---

## 3. Mathematical formulation

For a sequence of **L** unambiguous bases (A, C, G, T only), let *N_C*, *N_G* be
the numbers of C and G, and *N_CpG* the number of positions *i* with base *i* = C
and base *i+1* = G.

| Metric | Definition | Interpretation | If undefined |
|---|---|---|---|
| Length | number of symbols in the region | island must be long enough for statistics and biology | below minimum → fails |
| GC content | 100 × (N_C + N_G) / L | GC-rich vs AT-rich | L = 0 → `NA`, fails |
| Obs/Exp CpG | N_CpG × L / (N_C × N_G) | CpG retention relative to composition | N_C = 0 or N_G = 0 → `NA`, fails (never reported as 0 or ∞) |

**Why the formula estimates the expected CpG.** If bases were placed
independently with probabilities P(C) = N_C/L and P(G) = N_G/L, the chance that a
given adjacent pair is C→G is P(C)·P(G). With about L adjacent positions, the
expected count is L · (N_C/L) · (N_G/L) = N_C·N_G / L. Dividing the observed count
by this gives Obs/Exp = N_CpG · L / (N_C · N_G).

**Assumptions and limits.** (i) Independence of neighbouring bases; real DNA has
local correlations. (ii) Uniform composition; large regions are heterogeneous.
(iii) L adjacent pairs rather than L−1: a ≈ 0.5 % difference at 200 bp, negligible
here but not zero. (iv) Noisy for small counts. (v) Some published variants use
different but similar denominators, so results can differ slightly.

**Worked example.** `AACGCGTT`: L = 8, N_C = 2, N_G = 2, N_CpG = 2 (C3G4 and
C5G6); GC = 50 %; Obs/Exp = 2 × 8 / (2 × 2) = 4.0.

**Classification (defaults).** A region is a candidate if length ≥ 200 bp **and**
GC > 50 % **and** Obs/Exp > 0.6 **and** its ambiguous-base fraction ≤ 5 %. GC and
Obs/Exp use strict inequality; length is inclusive.

**Historical threshold vs universal rule.** The 200 bp / 50 % / 0.6 values are the
Gardiner-Garden & Frommer criteria, chosen to separate a particular set of human
sequences. They are not biological constants. Alternatives exist (Takai & Jones;
HMM-based and experiment-based definitions). All thresholds are adjustable on
the command line.

**Ambiguous bases.** IUPAC codes (N, R, Y, …) are accepted by default but never
converted to A/C/G/T. They are excluded from the GC and Obs/Exp denominators
(L counts only A/C/G/T), they cannot form part of a CpG (so `CNG` has no CpG),
and a window or region whose ambiguous fraction exceeds a limit (default 5 %) is
ineligible, because a handful of real bases among many N's would give meaningless
percentages. `--ambiguous reject` turns any ambiguity code into an error. Any
non-IUPAC symbol (digits, `U`, `-`, `*`) is always rejected with its position.

---

## 4. Detection algorithm and alternatives

### 4.1 Strategy 1: fixed-width sliding windows (implemented, with merging)
Score every window of width *W* starting every *S* bases; keep windows that pass
the GC and Obs/Exp criteria; merge overlapping or touching windows; re-measure
each merged region as a whole and keep it only if it still meets *all* criteria,
including minimum length.

* **Complexity.** Recounting each window is O(n/S · W). Using cumulative
  (prefix) counts of A, C, G, T and CpG makes every window O(1), so the scan is
  O(n + n/S) time and O(n) memory (five 32-bit integers per base).
  A 3 Mb random sequence scanned in about 2.6 s on the development machine.
* **Boundary effects.** A window passes if the island supplies enough of it, so
  qualifying windows hang over the true edges and merged boundaries overshoot,
  by up to about one window length, and are quantised to the step size.
* **Missed regions.** An island shorter than roughly half a window cannot lift a
  whole window above the thresholds and can be missed or diluted. A region whose
  best placement falls between window starts is partially seen.
* **Merging.** Windows that overlap or touch (next start ≤ current end) are
  merged. Windows separated by any gap stay separate. Re-checking the merged
  region guards against a union that fails although its parts pass.
* **End of sequence.** If the step grid does not land on the final base, one extra
  window flush with the 3' end is scanned so the tail is never skipped.
* **Why the defaults.** Window = 200 bp equals the minimum island length, so a
  single passing window is itself a valid-length island. Step = 10 bp gives a
  position error of ≈ 10 bp at modest cost. These are practical choices, not
  biological values.

### 4.2 Strategy 2: flexible region growing and merging (not implemented)
Seed at a small CpG-rich core, extend while criteria still hold, and trim ends
(Takai & Jones-style), or compute dynamic-programming/HMM segmentations.
* **Pros.** Boundaries are not tied to a grid and can follow the true edge more
  tightly; can be tuned to avoid Alu-like short repeats.
* **Cons.** More parameters (seed size, extension rule, trimming rule), greater
  implementation complexity, and order-dependent results unless carefully
  designed; harder to explain and test.

### 4.3 Choice
Strategy 1 is implemented because it is simple, deterministic, O(n), easy to
unit-test with hand-calculable cases, and its limitations are easy to state. Its
overshoot is documented rather than hidden. Boundary trimming is listed as a
future improvement.

### 4.4 Coordinates
Reports use **1-based, inclusive** coordinates. A Python slice `[s0, e0)`
(0-based, end-exclusive) becomes `(s0 + 1, e0)`. The reverse is
`sequence[start - 1 : end]`. The tests round-trip this conversion.

---

## 5. Software design

* **Single module, separated concerns:** validation/parsing → metrics →
  criteria → windows/merging → reporting/export → plotting → CLI. One file keeps
  the project easy to read and run; functions have single responsibilities.
* **Standard library only.** Plots are dependency-free SVG, so there is nothing
  to install and output is deterministic. A plotting library would add a
  dependency without changing the science.
* **Frozen dataclasses** for records, parameters, metrics and regions;
  parameters are validated on creation.
* **Deliberate errors:** `InvalidSequenceError`, `FastaFormatError`,
  `ParameterError` carry actionable messages (record name, line, position).
* **Determinism:** no randomness, timestamps or dict-order dependence in
  outputs; a test checks that two runs produce identical output.
* **unittest** was chosen over pytest because it ships with Python (no
  dependency) and is sufficient for this project.

---

## 6. Validation

### 6.1 Unit tests
77 `unittest` tests. Expected values were calculated by hand from the definitions
(each test carries its arithmetic in a comment), not by re-running the
implementation. Coverage includes: all-criteria pass; GC-rich with Obs/Exp = 0;
CpG-rich but too short; low GC (and exactly 50 % to test strict inequality);
lowercase; empty input; invalid characters (X, U, digits, gaps) with positions;
ambiguous N (exclusion from denominators; `CNG` is not a CpG; reject mode); zero
C or G; multi-record, multi-line, CRLF FASTA; malformed FASTA (sequence before
header, empty record, empty header, empty file, non-text file); overlapping CpGs
(`CGCG` = 2) and 1-based CpG positions; windows overlapping and merging; islands
at the very start/end of a sequence; sequence exactly one window long; sequence
shorter than the window; an end-of-sequence island found only through the
flush window; CSV/JSON export; SVG well-formedness; CLI exit codes.
A cross-check compares the O(1) prefix-count window metrics with direct counting
on a 600-bp pseudo-random sequence containing N, over many window placements.

### 6.2 Test of the tests (mutation check)
Four plausible bugs were injected one at a time (strict `>` changed to `>=` for
GC; end-flush window removed; Obs/Exp using total length instead of unambiguous
length; 1-based start off by one). The suite failed for each of them.

### 6.3 Synthetic validation (`data/sample.fasta`)
The sample sequences are **simulated** from base-composition models (realised on a
200-kb draw from each model: island-like GC ≈ 59 %, Obs/Exp ≈ 0.89; bulk-like
GC ≈ 36 %, Obs/Exp ≈ 0.28; decoy GC ≈ 61 %, Obs/Exp ≈ 0.29), with known island positions. The
generator (`scripts/generate_sample_fasta.py`) uses its own PRNG so the file is
reproducible across Python versions. Results with default settings:

| Record | Simulated truth | Predicted | Comment |
|---|---|---|---|
| synthetic_island_01 | 401–900 | 351–940 | overshoot 50 bp (5′) and 40 bp (3′) |
| synthetic_two_islands_02 | 301–650; 1151–1750 | 291–690; 1081–1840 | overshoot 10/40 bp and 70/90 bp |
| synthetic_gc_rich_cpg_depleted_03 | none | none | GC 60.9 %, Obs/Exp 0.222 → correctly rejected |
| synthetic_at_rich_background_04 | none | none | correct |
| synthetic_short_cpg_rich_05 | none (too short) | none | 150 bp; CpG-rich but fails length |
| synthetic_island_with_Ns_06 | 301–700 (N gap 451–480) | 241–460; 471–750 | island recovered but split at the N gap; flank overshoot 60/50 bp |

**Interpretation.** All islands are found, all negative controls are rejected, and
the decoy shows why Obs/Exp is needed. The overshoot is the boundary effect
described in §4.1. The N-gap split is a direct consequence of the rule that
windows containing > 5 % ambiguous bases are ineligible (30 N in a 200-bp window =
15 %). This is a design trade-off, not a biological finding.
**Caution:** this validates that the code does what it is designed to do on
simulated data. It does **not** measure accuracy on real genomes, because the
simulation was built from the same compositional ideas the detector uses.

### 6.4 What was and was not executed
The test suite, CLI sessions and 3 Mb timing were executed on Python 3.12.3.
Python 3.8–3.11 were **not** available, so the 3.8+ claim rests on avoiding newer
syntax and library features, not on execution. SVG files were checked for
well-formed XML but could not be rendered to an image in the development
environment, so their appearance should be checked in a browser.

---

## 7. Scientific honesty: statements the user must keep in mind
* Output = **computationally predicted candidate regions**, not confirmed CpG islands.
* Threshold-based detection is a useful educational and exploratory approach.
* Different definitions (200/50/0.6; 500/55/0.65; HMM; experimental) give different results.
* Fixed windows can miss or distort region boundaries.
* Sequence composition alone cannot establish promoter function, methylation
  status, gene expression or disease relevance; experimental evidence or external
  annotation is required.

## 8. Limitations and justified future improvements
| Limitation | Possible improvement and justification |
|---|---|
| Boundaries overshoot by up to ≈ a window | Trim region ends to the nearest CpG-dense core; reduces edge error shown in §6.3 |
| N gaps split regions | Option to bridge short gaps; trade-off between continuity and ambiguity |
| No repeat masking (Alu etc.) | Accept soft-masked (lowercase) input as a repeat filter; removes a known false-positive source |
| Both strands treated as one (CpG is palindromic, so counts are equal) | None needed for counts; stranded annotation matters only for gene context |
| Only 1987-style composition criteria | Add alternative presets (Takai-Jones) via a `--preset` flag |
| Memory is O(n) (≈ 20 bytes/base) | Chunked scanning for whole chromosomes |
| Single-threaded pure Python | Acceptable to tens of Mb; optional multiprocessing per record |
| No gene/annotation overlay | Report overlap with an annotation file (BED/GFF) to relate regions to promoters |
| No statistical significance of Obs/Exp | Permutation or Markov-model null for small windows |

## 9. References (to verify before citing)
* Gardiner-Garden M, Frommer M (1987). CpG islands in vertebrate genomes. *J Mol Biol* 196:261–282.
* Bird AP (1986). CpG-rich islands and the function of DNA methylation. *Nature* 321:209–213.
* Takai D, Jones PA (2002). Comprehensive analysis of CpG islands in human chromosomes 21 and 22. *PNAS* 99:3740–3745.
* Illingworth RS, Bird AP (2009). CpG islands: "a rough guide". *FEBS Lett* 583:1713–1720.
* Deaton AM, Bird A (2011). CpG islands and the regulation of transcription. *Genes Dev* 25:1010–1022.

(Bibliographic details were written from memory without internet access in the
development environment; check them against the originals before citing.)
