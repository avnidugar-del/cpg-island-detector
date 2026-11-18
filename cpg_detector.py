#!/usr/bin/env python3
"""CpG Island Detector: predict candidate CpG-rich regions in DNA sequences.

The program applies the historical Gardiner-Garden & Frommer (1987) style
criteria (length, GC content, observed/expected CpG ratio) to DNA sequences,
using a sliding-window scan followed by merging of qualifying windows.

IMPORTANT: outputs are *computationally predicted candidate regions* based on
sequence composition only. They are not experimentally confirmed CpG islands
and say nothing by themselves about promoter activity, methylation status,
gene expression or disease relevance.

Coordinate convention
---------------------
* Everything shown to the user (report, CSV, JSON, plots) uses 1-based,
  inclusive coordinates: a region "401-900" starts at the 401st base and ends
  at the 900th base, so its length is end - start + 1 = 500 bp.
* Internally, Python slicing is 0-based and end-exclusive. A region stored
  internally as the slice [start0, end0) is reported as
  (start0 + 1, end0). See :func:`to_one_based_inclusive`.

Requires Python 3.8 or newer and the standard library only.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from array import array
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple
from xml.sax.saxutils import escape as xml_escape

__version__ = "1.0.0"

# --------------------------------------------------------------------------
# Constants (no magic numbers elsewhere in the code)
# --------------------------------------------------------------------------

# Historical, commonly cited Gardiner-Garden & Frommer thresholds. These are
# conventions, not universal biological laws.
DEFAULT_MIN_LENGTH_BP = 200
DEFAULT_MIN_GC_PERCENT = 50.0
DEFAULT_MIN_OBS_EXP = 0.6

# Scan settings are practical choices, not biological constants.
DEFAULT_WINDOW_SIZE_BP = 200
DEFAULT_STEP_SIZE_BP = 10

# A window or region made mostly of N would give meaningless GC / Obs/Exp
# values from a handful of real bases, so heavily ambiguous windows are
# ineligible. 5 % is a conservative, adjustable default.
DEFAULT_MAX_AMBIGUOUS_FRACTION = 0.05

UNAMBIGUOUS_BASES = frozenset("ACGT")
# IUPAC ambiguity codes (N = any base, R = A/G, Y = C/T, ...).
IUPAC_AMBIGUOUS_BASES = frozenset("RYSWKMBDHVN")
ALL_ACCEPTED_SYMBOLS = UNAMBIGUOUS_BASES | IUPAC_AMBIGUOUS_BASES

MAX_INVALID_CHARS_REPORTED = 10
DIRECT_INPUT_RECORD_ID = "user_sequence"

TEXT_MAP_WIDTH = 60

# SVG plot layout
SVG_WIDTH = 900
SVG_HEIGHT = 540
SVG_MARGIN_LEFT = 70
SVG_MARGIN_RIGHT = 20
SVG_PANEL_HEIGHT = 150
SVG_MAX_PLOT_POINTS = 1200
SVG_OBS_EXP_AXIS_MAX = 1.5
SVG_RUG_BINS = 300

CSV_COLUMNS = (
    "sequence_id", "start_bp", "end_bp", "length_bp", "gc_percent",
    "obs_exp_cpg", "cpg_count", "c_count", "g_count", "a_count", "t_count",
    "ambiguous_count", "window_size_bp", "step_size_bp", "min_length_bp",
    "min_gc_percent", "min_obs_exp",
)

DISCLAIMER = (
    "Computationally predicted candidate regions based on sequence "
    "composition only; not experimentally confirmed CpG islands."
)


# --------------------------------------------------------------------------
# Exceptions
# --------------------------------------------------------------------------

class CpGDetectorError(Exception):
    """Base class for all deliberate, user-facing errors."""


class InvalidSequenceError(CpGDetectorError):
    """A sequence contains characters that are not valid DNA symbols."""


class FastaFormatError(CpGDetectorError):
    """A FASTA file or text is malformed."""


class ParameterError(CpGDetectorError, ValueError):
    """A detection parameter has an invalid value."""


# --------------------------------------------------------------------------
# Data structures
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SequenceRecord:
    """One named, validated, upper-case DNA sequence."""

    identifier: str
    description: str
    sequence: str


@dataclass(frozen=True)
class DetectionParameters:
    """Detection criteria and scan settings (validated on creation)."""

    min_length_bp: int = DEFAULT_MIN_LENGTH_BP
    min_gc_percent: float = DEFAULT_MIN_GC_PERCENT
    min_obs_exp: float = DEFAULT_MIN_OBS_EXP
    window_size_bp: int = DEFAULT_WINDOW_SIZE_BP
    step_size_bp: int = DEFAULT_STEP_SIZE_BP
    max_ambiguous_fraction: float = DEFAULT_MAX_AMBIGUOUS_FRACTION

    def __post_init__(self) -> None:
        if self.min_length_bp < 1:
            raise ParameterError("min-length must be at least 1 bp.")
        if not 0.0 <= self.min_gc_percent <= 100.0:
            raise ParameterError("min-gc must be between 0 and 100 (percent).")
        if self.min_obs_exp < 0.0:
            raise ParameterError("min-obs-exp must be 0 or greater.")
        if self.window_size_bp < 1:
            raise ParameterError("window-size must be at least 1 bp.")
        if self.step_size_bp < 1:
            raise ParameterError("step-size must be at least 1 bp.")
        if not 0.0 <= self.max_ambiguous_fraction <= 1.0:
            raise ParameterError(
                "max-ambiguous-fraction must be between 0 and 1.")


@dataclass(frozen=True)
class SequenceMetrics:
    """Composition metrics for one sequence or sub-sequence.

    ``length`` counts every symbol (including ambiguous ones). GC content and
    Obs/Exp use only the unambiguous A/C/G/T bases as their denominator.
    """

    length: int
    a_count: int
    c_count: int
    g_count: int
    t_count: int
    ambiguous_count: int
    cpg_count: int

    @property
    def unambiguous_length(self) -> int:
        return self.a_count + self.c_count + self.g_count + self.t_count

    @property
    def gc_percent(self) -> Optional[float]:
        """(G + C) / (A + C + G + T) * 100, or None if no A/C/G/T bases."""
        if self.unambiguous_length == 0:
            return None
        return 100.0 * (self.c_count + self.g_count) / self.unambiguous_length

    @property
    def obs_exp_cpg(self) -> Optional[float]:
        """CpG * L / (C * G), or None when C = 0 or G = 0 (undefined)."""
        if self.c_count == 0 or self.g_count == 0:
            return None
        return (self.cpg_count * self.unambiguous_length
                / (self.c_count * self.g_count))

    @property
    def ambiguous_fraction(self) -> float:
        return self.ambiguous_count / self.length if self.length else 0.0


@dataclass(frozen=True)
class CandidateRegion:
    """A predicted CpG-rich region; coordinates are 1-based and inclusive."""

    record_id: str
    start: int
    end: int
    metrics: SequenceMetrics

    @property
    def length(self) -> int:
        return self.end - self.start + 1


@dataclass(frozen=True)
class WindowPoint:
    """GC content and Obs/Exp of one scanned window (0-based, end-exclusive)."""

    start0: int
    end0: int
    gc_percent: Optional[float]
    obs_exp_cpg: Optional[float]


@dataclass
class RecordResult:
    """All analysis output for one sequence record."""

    record: SequenceRecord
    metrics: SequenceMetrics
    whole_sequence_failures: List[str]
    regions: List[CandidateRegion]
    profile: List[WindowPoint] = field(default_factory=list)


# --------------------------------------------------------------------------
# Coordinates
# --------------------------------------------------------------------------

def to_one_based_inclusive(start0: int, end0_exclusive: int) -> Tuple[int, int]:
    """Convert a Python slice [start0, end0) to 1-based inclusive (start, end)."""
    return start0 + 1, end0_exclusive


# --------------------------------------------------------------------------
# Input handling: validation and FASTA parsing
# --------------------------------------------------------------------------

def validate_sequence(raw_sequence: str, allow_ambiguous: bool = True,
                      label: str = "sequence") -> str:
    """Return the cleaned, upper-case sequence or raise InvalidSequenceError.

    * Whitespace (spaces, tabs, line breaks) is removed; reported positions
      refer to the sequence after whitespace removal.
    * Lower-case letters are accepted and converted to upper case (this does
      not change which nucleotide a symbol represents).
    * A, C, G, T are always valid. IUPAC ambiguity codes (N, R, Y, ...) are
      valid only if ``allow_ambiguous`` is True; they are never converted
      into A/C/G/T.
    * Anything else (digits, 'U', '-', '*', ...) is rejected, never repaired.
    """
    cleaned = "".join(raw_sequence.split()).upper()
    symbols_present = set(cleaned)

    invalid_symbols = symbols_present - ALL_ACCEPTED_SYMBOLS
    if invalid_symbols:
        found = [(index + 1, symbol) for index, symbol in enumerate(cleaned)
                 if symbol in invalid_symbols]
        shown = ", ".join(f"'{symbol}' at position {position}"
                          for position, symbol in
                          found[:MAX_INVALID_CHARS_REPORTED])
        more = ""
        if len(found) > MAX_INVALID_CHARS_REPORTED:
            more = (f" (showing first {MAX_INVALID_CHARS_REPORTED} of "
                    f"{len(found)})")
        hint = (" 'U' indicates RNA; convert to DNA first."
                if "U" in invalid_symbols else "")
        raise InvalidSequenceError(
            f"{label}: invalid character(s): {shown}{more}. Accepted "
            f"symbols are A, C, G, T (any case) and IUPAC ambiguity codes "
            f"{''.join(sorted(IUPAC_AMBIGUOUS_BASES))}.{hint}")

    if not allow_ambiguous:
        ambiguous_present = symbols_present & IUPAC_AMBIGUOUS_BASES
        if ambiguous_present:
            first_position = next(index + 1 for index, symbol
                                  in enumerate(cleaned)
                                  if symbol in ambiguous_present)
            raise InvalidSequenceError(
                f"{label}: ambiguous nucleotide code(s) "
                f"{''.join(sorted(ambiguous_present))} found (first at "
                f"position {first_position}) but ambiguous codes are "
                f"disallowed (--ambiguous reject).")
    return cleaned


def parse_fasta_lines(lines: Iterable[str], source: str = "<input>",
                      allow_ambiguous: bool = True) -> List[SequenceRecord]:
    """Parse FASTA text into validated records.

    Supports multiple records, multi-line sequences, blank lines, CRLF line
    endings and mixed case. The identifier is the first word of the header.
    """
    records: List[SequenceRecord] = []
    header: Optional[str] = None
    header_line_number = 0
    chunks: List[str] = []

    def finish_record() -> None:
        if header is None:
            return
        parts = header.split(None, 1)
        identifier = parts[0]
        description = parts[1] if len(parts) > 1 else ""
        if not chunks:
            raise FastaFormatError(
                f"{source}: record '{identifier}' (header on line "
                f"{header_line_number}) contains no sequence.")
        sequence = validate_sequence("".join(chunks), allow_ambiguous,
                                     label=f"{source}, record '{identifier}'")
        records.append(SequenceRecord(identifier, description, sequence))

    for line_number, raw_line in enumerate(lines, start=1):
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            finish_record()
            header = line[1:].strip()
            header_line_number = line_number
            chunks = []
            if not header:
                raise FastaFormatError(
                    f"{source}: empty FASTA header on line {line_number}; "
                    f"each record needs an identifier after '>'.")
        else:
            if header is None:
                raise FastaFormatError(
                    f"{source}: sequence data on line {line_number} appears "
                    f"before any '>' header line.")
            chunks.append(line)
    finish_record()

    if not records:
        raise FastaFormatError(f"{source}: no FASTA records found.")
    return records


def read_fasta_file(path: str, allow_ambiguous: bool = True
                    ) -> List[SequenceRecord]:
    """Read and parse a FASTA file (FileNotFoundError propagates)."""
    try:
        with open(path, "r", encoding="utf-8-sig") as handle:
            return parse_fasta_lines(handle, source=str(path),
                                     allow_ambiguous=allow_ambiguous)
    except UnicodeDecodeError as exc:
        raise FastaFormatError(
            f"{path}: not readable as UTF-8 text ({exc.reason}); "
            f"is this really a FASTA file?") from exc


# --------------------------------------------------------------------------
# Sequence analysis: metrics and criteria
# --------------------------------------------------------------------------

def compute_metrics(sequence: str, allow_ambiguous: bool = True
                    ) -> SequenceMetrics:
    """Compute composition metrics for a DNA sequence.

    CpG dinucleotides are counted as adjacent C->G pairs on the given strand.
    Because "CG" cannot overlap itself, ``str.count("CG")`` counts every
    occurrence (e.g. "CGCG" contains 2). Note that "GC" is not a CpG, and a
    C and G separated by an ambiguous base (e.g. "CNG") do not form a CpG.
    """
    cleaned = validate_sequence(sequence, allow_ambiguous)
    a_count = cleaned.count("A")
    c_count = cleaned.count("C")
    g_count = cleaned.count("G")
    t_count = cleaned.count("T")
    return SequenceMetrics(
        length=len(cleaned),
        a_count=a_count, c_count=c_count, g_count=g_count, t_count=t_count,
        ambiguous_count=len(cleaned) - (a_count + c_count + g_count + t_count),
        cpg_count=cleaned.count("CG"))


def find_cpg_positions(sequence: str) -> List[int]:
    """Return the 1-based position of the C of every CpG dinucleotide."""
    cleaned = validate_sequence(sequence)
    positions: List[int] = []
    index = cleaned.find("CG")
    while index != -1:
        positions.append(index + 1)
        index = cleaned.find("CG", index + 1)
    return positions


def evaluate_criteria(metrics: SequenceMetrics, params: DetectionParameters,
                      check_length: bool = True) -> List[str]:
    """Return the reasons the metrics fail the criteria (empty = passes).

    GC content and Obs/Exp must be *strictly greater* than their thresholds;
    length must be at least the minimum. Undefined values (zero denominator
    or no A/C/G/T bases) never pass.
    """
    failures: List[str] = []
    if check_length and metrics.length < params.min_length_bp:
        failures.append(f"length {metrics.length} bp is below the minimum "
                        f"{params.min_length_bp} bp")

    gc_percent = metrics.gc_percent
    if gc_percent is None:
        failures.append("GC content is undefined (no A/C/G/T bases)")
    elif not gc_percent > params.min_gc_percent:
        failures.append(f"GC content {gc_percent:.2f}% is not above "
                        f"{params.min_gc_percent:g}%")

    obs_exp = metrics.obs_exp_cpg
    if obs_exp is None:
        failures.append("Obs/Exp CpG is undefined (zero C or zero G count)")
    elif not obs_exp > params.min_obs_exp:
        failures.append(f"Obs/Exp CpG {obs_exp:.3f} is not above "
                        f"{params.min_obs_exp:g}")

    if metrics.ambiguous_fraction > params.max_ambiguous_fraction:
        failures.append(f"ambiguous-base fraction "
                        f"{metrics.ambiguous_fraction:.3f} exceeds "
                        f"{params.max_ambiguous_fraction:g}")
    return failures


# --------------------------------------------------------------------------
# Region detection: sliding windows + merging
# --------------------------------------------------------------------------

class PrefixCounts:
    """Cumulative base and CpG counts giving O(1) metrics for any window.

    ``cpg[k]`` is the number of CpG dinucleotides whose C lies at 0-based
    index < k. The CpGs fully inside the slice [start, end) have their C at
    indices start .. end-2, so their number is ``cpg[end - 1] - cpg[start]``.
    """

    def __init__(self, sequence: str) -> None:
        length = len(sequence)
        self.sequence_length = length
        self._a = array("i", [0]) * (length + 1)
        self._c = array("i", [0]) * (length + 1)
        self._g = array("i", [0]) * (length + 1)
        self._t = array("i", [0]) * (length + 1)
        self._cpg = array("i", [0]) * (length + 1)
        a_total = c_total = g_total = t_total = cpg_total = 0
        for index, base in enumerate(sequence):
            if base == "A":
                a_total += 1
            elif base == "C":
                c_total += 1
                # This C starts a CpG if the next base is G.
                if index + 1 < length and sequence[index + 1] == "G":
                    cpg_total += 1
            elif base == "G":
                g_total += 1
            elif base == "T":
                t_total += 1
            self._a[index + 1] = a_total
            self._c[index + 1] = c_total
            self._g[index + 1] = g_total
            self._t[index + 1] = t_total
            self._cpg[index + 1] = cpg_total

    def metrics(self, start: int, end: int) -> SequenceMetrics:
        """Metrics for the 0-based, end-exclusive slice [start, end)."""
        a_count = self._a[end] - self._a[start]
        c_count = self._c[end] - self._c[start]
        g_count = self._g[end] - self._g[start]
        t_count = self._t[end] - self._t[start]
        window_length = end - start
        cpg_count = (self._cpg[end - 1] - self._cpg[start]
                     if window_length >= 2 else 0)
        return SequenceMetrics(
            length=window_length, a_count=a_count, c_count=c_count,
            g_count=g_count, t_count=t_count,
            ambiguous_count=window_length - (a_count + c_count + g_count
                                             + t_count),
            cpg_count=cpg_count)


def generate_window_bounds(sequence_length: int, params: DetectionParameters
                           ) -> Iterator[Tuple[int, int]]:
    """Yield 0-based, end-exclusive window bounds covering the sequence.

    * Sequences shorter than ``min_length_bp`` cannot contain a qualifying
      region, so no windows are produced.
    * A sequence no longer than the window is scanned as one window.
    * Otherwise windows start at 0, step, 2*step, ...; one extra window
      flush with the sequence end is added when the step grid would not land
      on it, so the 3' end is never left unscanned.
    """
    if sequence_length < params.min_length_bp:
        return
    window = params.window_size_bp
    if sequence_length <= window:
        yield 0, sequence_length
        return
    last_start = sequence_length - window
    start = 0
    while start <= last_start:
        yield start, start + window
        start += params.step_size_bp
    if start - params.step_size_bp != last_start:
        yield last_start, sequence_length


def merge_intervals(intervals: Iterable[Tuple[int, int]]
                    ) -> List[Tuple[int, int]]:
    """Merge overlapping or touching [start, end) intervals."""
    merged: List[Tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def analyse_record(record: SequenceRecord, params: DetectionParameters,
                   include_profile: bool = False) -> RecordResult:
    """Scan one record and return whole-sequence metrics and candidate regions.

    Steps: (1) score every window with the criteria (length excluded, since a
    window may be shorter than the minimum island length); (2) merge
    overlapping or touching qualifying windows; (3) re-measure each merged
    region as a whole and keep it only if the region itself meets *all*
    criteria, including minimum length. Step 3 matters because the union of
    qualifying windows does not automatically qualify.
    """
    sequence = record.sequence
    whole_metrics = compute_metrics(sequence)
    whole_failures = evaluate_criteria(whole_metrics, params)

    prefix = PrefixCounts(sequence)
    qualifying: List[Tuple[int, int]] = []
    profile: List[WindowPoint] = []
    for start0, end0 in generate_window_bounds(len(sequence), params):
        window_metrics = prefix.metrics(start0, end0)
        if include_profile:
            profile.append(WindowPoint(start0, end0, window_metrics.gc_percent,
                                       window_metrics.obs_exp_cpg))
        if not evaluate_criteria(window_metrics, params, check_length=False):
            qualifying.append((start0, end0))

    regions: List[CandidateRegion] = []
    for start0, end0 in merge_intervals(qualifying):
        region_metrics = prefix.metrics(start0, end0)
        if evaluate_criteria(region_metrics, params):
            continue
        start, end = to_one_based_inclusive(start0, end0)
        regions.append(CandidateRegion(record.identifier, start, end,
                                       region_metrics))
    return RecordResult(record, whole_metrics, whole_failures, regions,
                        profile)


def detect_regions(record_id: str, sequence: str, params: DetectionParameters
                   ) -> List[CandidateRegion]:
    """Convenience wrapper: validate a raw sequence and return its regions."""
    record = SequenceRecord(record_id, "", validate_sequence(sequence))
    return analyse_record(record, params).regions


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def _format_number(value: Optional[float], decimals: int) -> str:
    return "NA" if value is None else f"{value:.{decimals}f}"


def render_text_map(sequence_length: int, regions: Sequence[CandidateRegion],
                    width: int = TEXT_MAP_WIDTH) -> str:
    """One-line map: '#' where a column overlaps a region, '.' elsewhere."""
    if sequence_length == 0:
        return ""
    width = min(width, sequence_length)
    cells = []
    for column in range(width):
        low = column * sequence_length // width
        high = max(low + 1, (column + 1) * sequence_length // width)
        hit = any(region.start - 1 < high and region.end > low
                  for region in regions)
        cells.append("#" if hit else ".")
    return "".join(cells)


def format_report(results: Sequence[RecordResult],
                  params: DetectionParameters) -> str:
    """Build the human-readable report (deterministic; no timestamps)."""
    lines = [
        f"CpG Island Detector v{__version__}",
        DISCLAIMER,
        (f"Criteria: length >= {params.min_length_bp} bp, "
         f"GC > {params.min_gc_percent:g} %, "
         f"Obs/Exp CpG > {params.min_obs_exp:g}, "
         f"ambiguous fraction <= {params.max_ambiguous_fraction:g}"),
        (f"Scan: window {params.window_size_bp} bp, "
         f"step {params.step_size_bp} bp. "
         f"Coordinates: 1-based, inclusive."),
        "",
    ]
    total_regions = 0
    for result in results:
        metrics = result.metrics
        total_regions += len(result.regions)
        lines.append(f"=== {result.record.identifier} "
                     f"({metrics.length} bp) ===")
        lines.append(
            f"Composition: A={metrics.a_count} C={metrics.c_count} "
            f"G={metrics.g_count} T={metrics.t_count} "
            f"ambiguous={metrics.ambiguous_count}")
        lines.append(
            f"Whole sequence: GC {_format_number(metrics.gc_percent, 2)} %, "
            f"CpG {metrics.cpg_count}, "
            f"Obs/Exp {_format_number(metrics.obs_exp_cpg, 3)}")
        if result.whole_sequence_failures:
            lines.append("Whole sequence meets all criteria: no")
            for reason in result.whole_sequence_failures:
                lines.append(f"  - {reason}")
        else:
            lines.append("Whole sequence meets all criteria: yes")
        lines.append(f"Candidate regions: {len(result.regions)}")
        if result.regions:
            lines.append(f"  {'#':>2} {'Start':>9} {'End':>9} {'Length':>7} "
                         f"{'GC%':>7} {'Obs/Exp':>8} {'CpG':>5}")
            for number, region in enumerate(result.regions, start=1):
                region_metrics = region.metrics
                lines.append(
                    f"  {number:>2} {region.start:>9} {region.end:>9} "
                    f"{region.length:>7} "
                    f"{_format_number(region_metrics.gc_percent, 2):>7} "
                    f"{_format_number(region_metrics.obs_exp_cpg, 3):>8} "
                    f"{region_metrics.cpg_count:>5}")
        if metrics.length:
            text_map = render_text_map(metrics.length, result.regions)
            bases_per_cell = metrics.length / len(text_map)
            lines.append(f"Map (1 char = {bases_per_cell:.1f} bp, "
                         f"# = candidate region):")
            lines.append(f"  |{text_map}|")
        lines.append("")
    lines.append(f"Total: {len(results)} sequence(s), "
                 f"{total_regions} candidate region(s).")
    return "\n".join(lines)


def _region_row(region: CandidateRegion, params: DetectionParameters
                ) -> Dict[str, object]:
    metrics = region.metrics
    return {
        "sequence_id": region.record_id,
        "start_bp": region.start,
        "end_bp": region.end,
        "length_bp": region.length,
        "gc_percent": _format_number(metrics.gc_percent, 2),
        "obs_exp_cpg": _format_number(metrics.obs_exp_cpg, 3),
        "cpg_count": metrics.cpg_count,
        "c_count": metrics.c_count,
        "g_count": metrics.g_count,
        "a_count": metrics.a_count,
        "t_count": metrics.t_count,
        "ambiguous_count": metrics.ambiguous_count,
        "window_size_bp": params.window_size_bp,
        "step_size_bp": params.step_size_bp,
        "min_length_bp": params.min_length_bp,
        "min_gc_percent": params.min_gc_percent,
        "min_obs_exp": params.min_obs_exp,
    }


def write_csv(results: Sequence[RecordResult], params: DetectionParameters,
              path: str) -> None:
    """Write one CSV row per candidate region (header only if none found)."""
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS,
                                lineterminator="\n")
        writer.writeheader()
        for result in results:
            for region in result.regions:
                writer.writerow(_region_row(region, params))


def write_json(results: Sequence[RecordResult], params: DetectionParameters,
               path: str) -> None:
    """Write parameters, per-sequence summaries and regions as JSON."""
    document = {
        "tool": "cpg_detector",
        "version": __version__,
        "disclaimer": DISCLAIMER,
        "coordinate_system": "1-based, inclusive",
        "parameters": {
            "min_length_bp": params.min_length_bp,
            "min_gc_percent": params.min_gc_percent,
            "min_obs_exp": params.min_obs_exp,
            "window_size_bp": params.window_size_bp,
            "step_size_bp": params.step_size_bp,
            "max_ambiguous_fraction": params.max_ambiguous_fraction,
        },
        "records": [],
    }
    for result in results:
        metrics = result.metrics
        document["records"].append({
            "sequence_id": result.record.identifier,
            "length_bp": metrics.length,
            "gc_percent": metrics.gc_percent,
            "obs_exp_cpg": metrics.obs_exp_cpg,
            "cpg_count": metrics.cpg_count,
            "ambiguous_count": metrics.ambiguous_count,
            "whole_sequence_meets_criteria":
                not result.whole_sequence_failures,
            "regions": [{
                "start_bp": region.start,
                "end_bp": region.end,
                "length_bp": region.length,
                "gc_percent": region.metrics.gc_percent,
                "obs_exp_cpg": region.metrics.obs_exp_cpg,
                "cpg_count": region.metrics.cpg_count,
            } for region in result.regions],
        })
    output_path = Path(path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as handle:
        json.dump(document, handle, indent=2)
        handle.write("\n")


def export_results(results: Sequence[RecordResult],
                   params: DetectionParameters, path: str) -> None:
    """Export to JSON if the path ends in .json, otherwise to CSV."""
    if path.lower().endswith(".json"):
        write_json(results, params, path)
    else:
        write_csv(results, params, path)


# --------------------------------------------------------------------------
# Visualisation (dependency-free SVG)
# --------------------------------------------------------------------------

def _polyline_markup(points: Sequence[Optional[Tuple[float, float]]],
                     colour: str) -> List[str]:
    """SVG polylines for runs of defined points; None breaks the line."""
    elements: List[str] = []
    run: List[Tuple[float, float]] = []
    for point in list(points) + [None]:
        if point is not None:
            run.append(point)
            continue
        if run:
            if len(run) == 1:
                run = run * 2  # a round-capped zero-length line draws a dot
            coordinates = " ".join(f"{x:.1f},{y:.1f}" for x, y in run)
            elements.append(
                f'<polyline fill="none" stroke="{colour}" stroke-width="1.5" '
                f'stroke-linecap="round" stroke-linejoin="round" '
                f'points="{coordinates}"/>')
            run = []
    return elements


def build_svg_profile(result: RecordResult, params: DetectionParameters) -> str:
    """Return SVG markup: GC%, Obs/Exp, predicted regions and CpG density."""
    sequence = result.record.sequence
    length = len(sequence)
    plot_width = SVG_WIDTH - SVG_MARGIN_LEFT - SVG_MARGIN_RIGHT
    gc_top, oe_top = 50, 240
    region_top, rug_top, axis_y = 430, 460, 500

    def x_of(position0: float) -> float:
        return SVG_MARGIN_LEFT + plot_width * position0 / length

    def gc_y(value: float) -> float:
        return gc_top + SVG_PANEL_HEIGHT * (1 - value / 100.0)

    def oe_y(value: float) -> float:
        capped = min(value, SVG_OBS_EXP_AXIS_MAX)
        return oe_top + SVG_PANEL_HEIGHT * (1 - capped / SVG_OBS_EXP_AXIS_MAX)

    title = xml_escape(f"{result.record.identifier} ({length} bp)")
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{SVG_WIDTH}" '
        f'height="{SVG_HEIGHT}" viewBox="0 0 {SVG_WIDTH} {SVG_HEIGHT}" '
        f'font-family="sans-serif" font-size="12">',
        f'<title>CpG profile: {title}</title>',
        f'<rect width="{SVG_WIDTH}" height="{SVG_HEIGHT}" fill="white"/>',
        f'<text x="{SVG_MARGIN_LEFT}" y="22" font-size="14" '
        f'font-weight="bold">{title}</text>',
        f'<text x="{SVG_MARGIN_LEFT}" y="38" fill="#555">'
        f'Predicted candidate regions (window {params.window_size_bp} bp, '
        f'step {params.step_size_bp} bp)</text>',
    ]

    # Shaded candidate regions behind both panels.
    for region in result.regions:
        x_start = x_of(region.start - 1)
        x_width = max(1.0, x_of(region.end) - x_start)
        parts.append(
            f'<rect x="{x_start:.1f}" y="{gc_top}" width="{x_width:.1f}" '
            f'height="{oe_top + SVG_PANEL_HEIGHT - gc_top}" '
            f'fill="#cfe3f7" opacity="0.7"/>')

    # Panel frames, labels and threshold lines.
    for top, label in ((gc_top, "GC content (%)"),
                       (oe_top, "Obs/Exp CpG")):
        parts.append(f'<rect x="{SVG_MARGIN_LEFT}" y="{top}" '
                     f'width="{plot_width}" height="{SVG_PANEL_HEIGHT}" '
                     f'fill="none" stroke="#999"/>')
        parts.append(f'<text x="{SVG_MARGIN_LEFT}" y="{top - 6}">'
                     f'{label}</text>')
    gc_threshold_y = gc_y(params.min_gc_percent)
    oe_threshold_y = oe_y(params.min_obs_exp)
    for y_value, text in ((gc_threshold_y, f"{params.min_gc_percent:g}"),
                          (oe_threshold_y, f"{params.min_obs_exp:g}")):
        parts.append(
            f'<line x1="{SVG_MARGIN_LEFT}" x2="{SVG_WIDTH - SVG_MARGIN_RIGHT}"'
            f' y1="{y_value:.1f}" y2="{y_value:.1f}" stroke="#c0392b" '
            f'stroke-dasharray="5,4"/>')
        parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" y="{y_value + 4:.1f}" '
                     f'text-anchor="end" fill="#c0392b">{text}</text>')
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" y="{gc_top + 4}" '
                 f'text-anchor="end">100</text>')
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" '
                 f'y="{gc_top + SVG_PANEL_HEIGHT + 4}" '
                 f'text-anchor="end">0</text>')
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" y="{oe_top + 4}" '
                 f'text-anchor="end">{SVG_OBS_EXP_AXIS_MAX:g}+</text>')
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" '
                 f'y="{oe_top + SVG_PANEL_HEIGHT + 4}" '
                 f'text-anchor="end">0</text>')

    # Window profile curves (downsampled for long sequences).
    stride = max(1, len(result.profile) // SVG_MAX_PLOT_POINTS)
    sampled = result.profile[::stride]
    gc_points: List[Optional[Tuple[float, float]]] = []
    oe_points: List[Optional[Tuple[float, float]]] = []
    for point in sampled:
        centre_x = x_of((point.start0 + point.end0) / 2)
        gc_points.append(None if point.gc_percent is None
                         else (centre_x, gc_y(point.gc_percent)))
        oe_points.append(None if point.obs_exp_cpg is None
                         else (centre_x, oe_y(point.obs_exp_cpg)))
    parts.extend(_polyline_markup(gc_points, "#1f6fb2"))
    parts.extend(_polyline_markup(oe_points, "#2e8b57"))

    # Region track.
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" y="{region_top + 13}" '
                 f'text-anchor="end">Regions</text>')
    parts.append(f'<rect x="{SVG_MARGIN_LEFT}" y="{region_top}" '
                 f'width="{plot_width}" height="18" fill="#f4f4f4" '
                 f'stroke="#ccc"/>')
    for region in result.regions:
        x_start = x_of(region.start - 1)
        x_width = max(1.0, x_of(region.end) - x_start)
        parts.append(
            f'<rect x="{x_start:.1f}" y="{region_top}" '
            f'width="{x_width:.1f}" height="18" fill="#2a6fdb">'
            f'<title>{region.start}-{region.end} '
            f'({region.length} bp)</title></rect>')

    # CpG density strip: opacity proportional to CpGs per bin.
    parts.append(f'<text x="{SVG_MARGIN_LEFT - 6}" y="{rug_top + 13}" '
                 f'text-anchor="end">CpG density</text>')
    bin_count = min(SVG_RUG_BINS, length)
    counts = [0] * bin_count
    for position in find_cpg_positions(sequence):
        counts[min(bin_count - 1, (position - 1) * bin_count // length)] += 1
    busiest = max(counts) if counts else 0
    parts.append(f'<rect x="{SVG_MARGIN_LEFT}" y="{rug_top}" '
                 f'width="{plot_width}" height="18" fill="none" '
                 f'stroke="#ccc"/>')
    if busiest:
        bin_width = plot_width / bin_count
        for index, count in enumerate(counts):
            if count:
                parts.append(
                    f'<rect x="{SVG_MARGIN_LEFT + index * bin_width:.2f}" '
                    f'y="{rug_top}" width="{bin_width + 0.3:.2f}" '
                    f'height="18" fill="#7b2d8e" '
                    f'opacity="{count / busiest:.2f}"/>')

    # Position axis.
    parts.append(f'<line x1="{SVG_MARGIN_LEFT}" x2="{SVG_WIDTH - SVG_MARGIN_RIGHT}"'
                 f' y1="{axis_y}" y2="{axis_y}" stroke="#333"/>')
    tick_count = 5
    for tick in range(tick_count + 1):
        position = round(length * tick / tick_count)
        x_value = x_of(position)
        parts.append(f'<line x1="{x_value:.1f}" x2="{x_value:.1f}" '
                     f'y1="{axis_y}" y2="{axis_y + 5}" stroke="#333"/>')
        parts.append(f'<text x="{x_value:.1f}" y="{axis_y + 19}" '
                     f'text-anchor="middle">{position}</text>')
    parts.append(f'<text x="{SVG_MARGIN_LEFT + plot_width / 2:.0f}" '
                 f'y="{axis_y + 36}" text-anchor="middle">'
                 f'Position (bp)</text>')
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def write_svg_profile(result: RecordResult, params: DetectionParameters,
                      path: Path) -> None:
    """Write the SVG profile for one record to ``path``."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_svg_profile(result, params), encoding="utf-8")


def safe_file_stem(identifier: str) -> str:
    """Make a record identifier safe to use inside a file name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", identifier).strip("._") or "record"


# --------------------------------------------------------------------------
# Command-line interface
# --------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cpg_detector.py",
        description=("Predict candidate CpG-rich regions in DNA sequences "
                     "using sequence-composition criteria. " + DISCLAIMER),
        epilog=("Coordinates are 1-based and inclusive. Defaults follow the "
                "historical Gardiner-Garden & Frommer criteria, which are "
                "conventions rather than universal biological rules."))
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("-s", "--sequence", metavar="DNA",
                        help="DNA sequence given directly on the command line")
    source.add_argument("-f", "--fasta", metavar="FILE",
                        help="FASTA file (one or more records)")
    source.add_argument("-i", "--interactive", action="store_true",
                        help="type or paste a DNA sequence at a prompt")

    criteria = parser.add_argument_group("detection criteria")
    criteria.add_argument("--min-length", type=int,
                          default=DEFAULT_MIN_LENGTH_BP, metavar="BP",
                          help="minimum region length (default: %(default)s)")
    criteria.add_argument("--min-gc", type=float,
                          default=DEFAULT_MIN_GC_PERCENT, metavar="PERCENT",
                          help="GC%% must be strictly above this "
                               "(default: %(default)s)")
    criteria.add_argument("--min-obs-exp", type=float,
                          default=DEFAULT_MIN_OBS_EXP, metavar="RATIO",
                          help="Obs/Exp CpG must be strictly above this "
                               "(default: %(default)s)")
    criteria.add_argument("--max-ambiguous-fraction", type=float,
                          default=DEFAULT_MAX_AMBIGUOUS_FRACTION,
                          metavar="FRACTION",
                          help="largest allowed fraction of ambiguous bases "
                               "in a window/region (default: %(default)s)")

    scan = parser.add_argument_group("scan settings")
    scan.add_argument("--window-size", type=int,
                      default=DEFAULT_WINDOW_SIZE_BP, metavar="BP",
                      help="sliding-window width (default: %(default)s)")
    scan.add_argument("--step-size", type=int, default=DEFAULT_STEP_SIZE_BP,
                      metavar="BP",
                      help="distance between window starts "
                           "(default: %(default)s)")

    parser.add_argument("--ambiguous", choices=("allow", "reject"),
                        default="allow",
                        help="allow IUPAC ambiguity codes such as N "
                             "(excluded from GC/ObsExp denominators) or "
                             "reject them (default: %(default)s)")
    parser.add_argument("-o", "--output", metavar="FILE",
                        help="export regions; .json gives JSON, anything "
                             "else gives CSV")
    parser.add_argument("--plot-dir", metavar="DIR",
                        help="write one SVG profile per sequence to DIR")
    parser.add_argument("--version", action="version",
                        version=f"%(prog)s {__version__}")
    return parser


def load_records(args: argparse.Namespace, allow_ambiguous: bool
                 ) -> List[SequenceRecord]:
    """Obtain validated records from whichever input source was chosen."""
    if args.fasta is not None:
        return read_fasta_file(args.fasta, allow_ambiguous)
    if args.sequence is not None:
        raw_sequence = args.sequence
    else:
        try:
            raw_sequence = input("Enter DNA sequence: ")
        except EOFError as exc:
            raise CpGDetectorError("no sequence was entered.") from exc
    sequence = validate_sequence(raw_sequence, allow_ambiguous,
                                 label="input sequence")
    return [SequenceRecord(DIRECT_INPUT_RECORD_ID, "", sequence)]


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Program entry point; returns a process exit code."""
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        params = DetectionParameters(
            min_length_bp=args.min_length, min_gc_percent=args.min_gc,
            min_obs_exp=args.min_obs_exp, window_size_bp=args.window_size,
            step_size_bp=args.step_size,
            max_ambiguous_fraction=args.max_ambiguous_fraction)
    except ParameterError as exc:
        parser.error(str(exc))  # prints usage + message, exits with code 2

    try:
        records = load_records(args, allow_ambiguous=args.ambiguous == "allow")
        results = [analyse_record(record, params,
                                  include_profile=args.plot_dir is not None)
                   for record in records]
        print(format_report(results, params))
        if args.output:
            export_results(results, params, args.output)
            print(f"Results written to {args.output}")
        if args.plot_dir:
            for index, result in enumerate(results, start=1):
                plot_path = Path(args.plot_dir) / (
                    f"{index:03d}_{safe_file_stem(result.record.identifier)}"
                    f".svg")
                write_svg_profile(result, params, plot_path)
                print(f"Plot written to {plot_path}")
    except FileNotFoundError as exc:
        print(f"error: file not found: {exc.filename}", file=sys.stderr)
        return 1
    except CpGDetectorError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    except OSError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
