"""Unit tests for cpg_detector.

Expected values were calculated by hand from the definitions (see the comment
beside each case), not by calling the implementation. Run from the project
root with:  python -m unittest discover -s tests -v
"""

import contextlib
import csv
import io
import json
import sys
import tempfile
import unittest
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))

import cpg_detector as cpg  # noqa: E402

DEFAULTS = cpg.DetectionParameters()

# Hand-designed building blocks (all lengths and counts verified by hand).
# ISLAND_UNIT = "CGCGCGAT": C=3 G=3 A=1 T=1, CpG=3 (the unit junction T|C
# creates no extra CpG). GC = 6/8 = 75 %.
ISLAND_UNIT = "CGCGCGAT"
# FLANK_UNIT = "AT": no C or G at all.
FLANK_UNIT = "AT"


def island(length_bp: int) -> str:
    assert length_bp % len(ISLAND_UNIT) == 0
    return ISLAND_UNIT * (length_bp // len(ISLAND_UNIT))


def flank(length_bp: int) -> str:
    assert length_bp % len(FLANK_UNIT) == 0
    return FLANK_UNIT * (length_bp // len(FLANK_UNIT))


def lcg_sequence(length: int, seed: int = 7, alphabet: str = "ACGTN") -> str:
    """Tiny deterministic pseudo-random sequence (for property checks)."""
    state, symbols = seed, []
    for _ in range(length):
        state = (state * 1103515245 + 12345) % 2 ** 31
        symbols.append(alphabet[(state >> 16) % len(alphabet)])
    return "".join(symbols)


class TestValidation(unittest.TestCase):
    def test_lowercase_is_accepted_and_uppercased(self):
        self.assertEqual(cpg.validate_sequence("acgt"), "ACGT")

    def test_whitespace_and_line_breaks_removed(self):
        self.assertEqual(cpg.validate_sequence("AC\nG T\r\n"), "ACGT")

    def test_invalid_character_reported_with_position(self):
        with self.assertRaises(cpg.InvalidSequenceError) as context:
            cpg.validate_sequence("ACGTX")
        message = str(context.exception)
        self.assertIn("'X' at position 5", message)

    def test_rna_u_is_rejected_not_converted(self):
        with self.assertRaises(cpg.InvalidSequenceError) as context:
            cpg.validate_sequence("ACGU")
        self.assertIn("'U' at position 4", str(context.exception))
        self.assertIn("RNA", str(context.exception))

    def test_digits_and_gaps_rejected(self):
        for bad in ("AC1G", "AC-G", "AC*G"):
            with self.assertRaises(cpg.InvalidSequenceError):
                cpg.validate_sequence(bad)

    def test_many_invalid_characters_are_truncated_in_message(self):
        with self.assertRaises(cpg.InvalidSequenceError) as context:
            cpg.validate_sequence("A" + "X" * 25)
        self.assertIn("first 10 of 25", str(context.exception))

    def test_ambiguous_allowed_by_default(self):
        self.assertEqual(cpg.validate_sequence("ACGNRY"), "ACGNRY")

    def test_ambiguous_can_be_rejected(self):
        with self.assertRaises(cpg.InvalidSequenceError) as context:
            cpg.validate_sequence("ACGNT", allow_ambiguous=False)
        self.assertIn("position 4", str(context.exception))


class TestMetrics(unittest.TestCase):
    def test_counts_of_known_sequence(self):
        # "AACGCGTT": A=2 C=2 G=2 T=2, CpG at (3,4) and (5,6) -> 2.
        metrics = cpg.compute_metrics("AACGCGTT")
        self.assertEqual((metrics.length, metrics.a_count, metrics.c_count,
                          metrics.g_count, metrics.t_count), (8, 2, 2, 2, 2))
        self.assertEqual(metrics.cpg_count, 2)
        self.assertEqual(metrics.gc_percent, 50.0)
        # Obs/Exp = 2 * 8 / (2 * 2) = 4.0
        self.assertEqual(metrics.obs_exp_cpg, 4.0)

    def test_lowercase_gives_same_metrics_as_uppercase(self):
        self.assertEqual(cpg.compute_metrics("cgcgcgat" * 25),
                         cpg.compute_metrics("CGCGCGAT" * 25))

    def test_overlapping_cpg_dinucleotides_are_all_counted(self):
        self.assertEqual(cpg.compute_metrics("CGCGCG").cpg_count, 3)
        self.assertEqual(cpg.compute_metrics("CGCG").cpg_count, 2)

    def test_gc_is_not_a_cpg_and_pairs_are_strand_specific(self):
        self.assertEqual(cpg.compute_metrics("GCGC").cpg_count, 1)  # C(2)G(3)
        self.assertEqual(cpg.compute_metrics("GGCC").cpg_count, 0)
        self.assertEqual(cpg.compute_metrics("CCGG").cpg_count, 1)

    def test_cpg_positions_are_one_based_at_the_c(self):
        # A1 A2 C3 G4 C5 G6 T7 T8 -> C of each CpG at 3 and 5.
        self.assertEqual(cpg.find_cpg_positions("AACGCGTT"), [3, 5])
        self.assertEqual(cpg.find_cpg_positions("ATAT"), [])
        self.assertEqual(cpg.find_cpg_positions("cg"), [1])

    def test_empty_sequence_is_safe(self):
        metrics = cpg.compute_metrics("")
        self.assertEqual(metrics.length, 0)
        self.assertIsNone(metrics.gc_percent)
        self.assertIsNone(metrics.obs_exp_cpg)
        failures = cpg.evaluate_criteria(metrics, DEFAULTS)
        self.assertTrue(any("below the minimum" in f for f in failures))
        self.assertTrue(any("GC content is undefined" in f for f in failures))

    def test_zero_c_or_g_gives_undefined_obs_exp_not_zero_division(self):
        self.assertIsNone(cpg.compute_metrics("ATATATAT").obs_exp_cpg)  # C=G=0
        self.assertIsNone(cpg.compute_metrics("CCCC").obs_exp_cpg)      # G=0
        self.assertIsNone(cpg.compute_metrics("GGGG").obs_exp_cpg)      # C=0
        self.assertEqual(cpg.compute_metrics("ATATATAT").gc_percent, 0.0)
        self.assertEqual(cpg.compute_metrics("CCCC").gc_percent, 100.0)

    def test_undefined_obs_exp_never_passes(self):
        failures = cpg.evaluate_criteria(cpg.compute_metrics("C" * 300),
                                         DEFAULTS)
        self.assertTrue(any("Obs/Exp CpG is undefined" in f
                            for f in failures))

    def test_ambiguous_bases_excluded_from_denominators(self):
        # "ACGNNT": A=1 C=1 G=1 T=1 N=2, length 6, unambiguous L=4, CpG=1.
        metrics = cpg.compute_metrics("ACGNNT")
        self.assertEqual(metrics.length, 6)
        self.assertEqual(metrics.ambiguous_count, 2)
        self.assertEqual(metrics.unambiguous_length, 4)
        self.assertEqual(metrics.gc_percent, 50.0)       # 2 / 4
        self.assertEqual(metrics.obs_exp_cpg, 4.0)       # 1 * 4 / (1 * 1)

    def test_n_between_c_and_g_is_not_a_cpg(self):
        self.assertEqual(cpg.compute_metrics("CNG").cpg_count, 0)

    def test_compute_metrics_never_hides_invalid_characters(self):
        with self.assertRaises(cpg.InvalidSequenceError):
            cpg.compute_metrics("ACGTZ")


class TestCriteria(unittest.TestCase):
    def failures(self, sequence):
        return cpg.evaluate_criteria(cpg.compute_metrics(sequence), DEFAULTS)

    def test_sequence_meeting_all_criteria(self):
        # island(200): C=75 G=75 CpG=75 L=200 -> GC 75 %,
        # Obs/Exp = 75 * 200 / (75 * 75) = 2.6667.
        sequence = island(200)
        metrics = cpg.compute_metrics(sequence)
        self.assertAlmostEqual(metrics.gc_percent, 75.0)
        self.assertAlmostEqual(metrics.obs_exp_cpg, 8 / 3)
        self.assertEqual(self.failures(sequence), [])

    def test_gc_rich_but_low_obs_exp(self):
        # G*100 + C*100: GC = 100 %, only junction is "GC" -> CpG = 0,
        # Obs/Exp = 0. Fails ONLY the Obs/Exp criterion.
        failures = self.failures("G" * 100 + "C" * 100)
        self.assertEqual(len(failures), 1)
        self.assertIn("Obs/Exp CpG 0.000", failures[0])

    def test_cpg_abundant_but_too_short(self):
        # "CG" * 50: 100 bp, GC 100 %, Obs/Exp = 50*100/(50*50) = 2.0.
        failures = self.failures("CG" * 50)
        self.assertEqual(len(failures), 1)
        self.assertIn("length 100 bp", failures[0])

    def test_insufficient_gc_content(self):
        # "CGATAT" * 34: 204 bp, C=34 G=34, GC = 68/204 = 33.3 %,
        # CpG = 34, Obs/Exp = 34*204/(34*34) = 6.0 (high) -> fails GC only.
        sequence = "CGATAT" * 34
        self.assertEqual(cpg.compute_metrics(sequence).obs_exp_cpg, 6.0)
        failures = self.failures(sequence)
        self.assertEqual(len(failures), 1)
        self.assertIn("GC content 33.33%", failures[0])

    def test_gc_exactly_at_threshold_is_rejected_strictly_greater(self):
        # "CGAT" * 50: GC = exactly 50 % -> not "> 50".
        failures = self.failures("CGAT" * 50)
        self.assertEqual(len(failures), 1)
        self.assertIn("GC content 50.00%", failures[0])

    def test_length_exactly_at_minimum_is_accepted(self):
        self.assertEqual(self.failures(island(200)), [])
        self.assertTrue(any("length" in f for f in self.failures(island(192))))

    def test_heavily_ambiguous_window_is_ineligible(self):
        sequence = island(200) + "N" * 100          # 100 / 300 = 0.333 N
        failures = self.failures(sequence)
        self.assertTrue(any("ambiguous-base fraction" in f
                            for f in failures))


class TestParameters(unittest.TestCase):
    def test_invalid_parameters_are_rejected(self):
        bad_settings = (
            {"min_length_bp": 0}, {"min_gc_percent": -1},
            {"min_gc_percent": 100.5}, {"min_obs_exp": -0.1},
            {"window_size_bp": 0}, {"step_size_bp": 0},
            {"max_ambiguous_fraction": 1.5},
        )
        for settings in bad_settings:
            with self.subTest(settings=settings):
                with self.assertRaises(cpg.ParameterError):
                    cpg.DetectionParameters(**settings)


class TestFasta(unittest.TestCase):
    def test_multiple_records_multiline_lowercase_crlf(self):
        text = (">seq1 first sequence\r\nACGT\r\nacgt\r\n\r\n"
                ">seq2\r\nTTTT\r\n>seq3 third\r\nNNAC\r\n")
        records = cpg.parse_fasta_lines(io.StringIO(text))
        self.assertEqual([r.identifier for r in records],
                         ["seq1", "seq2", "seq3"])
        self.assertEqual([r.sequence for r in records],
                         ["ACGTACGT", "TTTT", "NNAC"])
        self.assertEqual(records[0].description, "first sequence")
        self.assertEqual(records[1].description, "")

    def test_sequence_before_header_is_an_error(self):
        with self.assertRaises(cpg.FastaFormatError) as context:
            cpg.parse_fasta_lines(["ACGT\n", ">seq\n", "ACGT\n"])
        self.assertIn("line 1", str(context.exception))

    def test_record_without_sequence_is_an_error(self):
        with self.assertRaises(cpg.FastaFormatError) as context:
            cpg.parse_fasta_lines([">empty\n", ">next\n", "ACGT\n"])
        self.assertIn("'empty'", str(context.exception))

    def test_empty_header_is_an_error(self):
        with self.assertRaises(cpg.FastaFormatError):
            cpg.parse_fasta_lines([">\n", "ACGT\n"])

    def test_empty_input_is_an_error(self):
        with self.assertRaises(cpg.FastaFormatError):
            cpg.parse_fasta_lines([])
        with self.assertRaises(cpg.FastaFormatError):
            cpg.parse_fasta_lines(["\n", "\n"])

    def test_invalid_character_error_names_the_record(self):
        with self.assertRaises(cpg.InvalidSequenceError) as context:
            cpg.parse_fasta_lines([">good\n", "ACGT\n", ">bad\n", "AC9T\n"])
        self.assertIn("record 'bad'", str(context.exception))

    def test_ambiguous_rejected_in_fasta_when_requested(self):
        with self.assertRaises(cpg.InvalidSequenceError):
            cpg.parse_fasta_lines([">s\n", "ACNT\n"], allow_ambiguous=False)

    def test_read_file_round_trip_and_missing_file(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "test.fasta"
            path.write_text(">a\nACGT\n>b\nGGCC\n", encoding="utf-8")
            records = cpg.read_fasta_file(str(path))
            self.assertEqual(len(records), 2)
            with self.assertRaises(FileNotFoundError):
                cpg.read_fasta_file(str(Path(folder) / "missing.fasta"))

    def test_non_text_file_gives_clear_error(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "binary.fasta"
            path.write_bytes(b"\xff\xfe\x00\x80\x81")
            with self.assertRaises(cpg.FastaFormatError):
                cpg.read_fasta_file(str(path))


class TestWindowsAndMerging(unittest.TestCase):
    @staticmethod
    def params(window, step, **kwargs):
        return cpg.DetectionParameters(window_size_bp=window,
                                       step_size_bp=step, **kwargs)

    def test_window_bounds_cover_end_when_grid_misses_it(self):
        # n=450, window 200, step 100 -> starts 0,100,200 then flush start 250.
        bounds = list(cpg.generate_window_bounds(450, self.params(200, 100)))
        self.assertEqual(bounds, [(0, 200), (100, 300), (200, 400),
                                  (250, 450)])

    def test_window_bounds_when_grid_lands_exactly_on_end(self):
        # n=400, window 200, step 100 -> 0,100,200 and 200+200 == 400: no extra.
        bounds = list(cpg.generate_window_bounds(400, self.params(200, 100)))
        self.assertEqual(bounds, [(0, 200), (100, 300), (200, 400)])

    def test_sequence_equal_to_or_shorter_than_window_is_one_window(self):
        self.assertEqual(list(cpg.generate_window_bounds(
            200, self.params(200, 10))), [(0, 200)])
        self.assertEqual(list(cpg.generate_window_bounds(
            250, self.params(300, 10))), [(0, 250)])

    def test_sequence_shorter_than_min_length_has_no_windows(self):
        self.assertEqual(list(cpg.generate_window_bounds(
            199, self.params(200, 10))), [])
        self.assertEqual(list(cpg.generate_window_bounds(
            0, self.params(200, 10))), [])

    def test_merge_overlapping_touching_and_separate_intervals(self):
        merged = cpg.merge_intervals([(300, 500), (0, 100), (50, 150),
                                      (150, 200), (400, 450)])
        # (0,100)+(50,150)+(150,200) touch/overlap -> (0,200);
        # (300,500) swallows (400,450).
        self.assertEqual(merged, [(0, 200), (300, 500)])

    def test_prefix_window_metrics_match_direct_counting(self):
        # Two independent implementations must agree on every window,
        # including windows containing N and windows straddling CpGs.
        sequence = lcg_sequence(600)
        prefix = cpg.PrefixCounts(sequence)
        for start in range(0, 560, 37):
            for offset in (1, 2, 50, 211):
                end = min(600, start + offset)
                self.assertEqual(prefix.metrics(start, end),
                                 cpg.compute_metrics(sequence[start:end]),
                                 msg=f"window [{start},{end})")

    def test_cpg_straddling_a_window_boundary_is_not_counted_in_either(self):
        # "ACGT": CpG occupies indices 1-2. Window [0,2) = "AC", [2,4) = "GT".
        prefix = cpg.PrefixCounts("ACGT")
        self.assertEqual(prefix.metrics(0, 2).cpg_count, 0)
        self.assertEqual(prefix.metrics(2, 4).cpg_count, 0)
        self.assertEqual(prefix.metrics(1, 3).cpg_count, 1)


class TestRegionDetection(unittest.TestCase):
    @staticmethod
    def detect(sequence, window=200, step=100, **kwargs):
        params = cpg.DetectionParameters(window_size_bp=window,
                                         step_size_bp=step, **kwargs)
        return cpg.detect_regions("test", sequence, params)

    def test_known_island_with_overlapping_windows_merges_to_one_region(self):
        # 400 AT + 400 island + 400 AT = 1200 bp, island = 401-800.
        # Window (step 100) analysis by hand:
        #   [300,500) has 100 island bases -> GC 37.5 % -> fails
        #   [400,600), [500,700), [600,800) are inside the island -> pass
        #   [700,900) has 100 island bases -> fails
        # Passing windows overlap, merge to [400,800) -> 401-800.
        regions = self.detect(flank(400) + island(400) + flank(400))
        self.assertEqual(len(regions), 1)
        region = regions[0]
        self.assertEqual((region.start, region.end, region.length),
                         (401, 800, 400))
        # Region metrics: C=150 G=150 CpG=150 L=400.
        self.assertEqual(region.metrics.cpg_count, 150)
        self.assertAlmostEqual(region.metrics.gc_percent, 75.0)
        self.assertAlmostEqual(region.metrics.obs_exp_cpg,
                               150 * 400 / (150 * 150))

    def test_two_separated_islands_give_two_regions(self):
        # island 401-800 and island 1401-1800 with 600 bp of AT between.
        sequence = (flank(400) + island(400) + flank(600) + island(400)
                    + flank(400))
        regions = self.detect(sequence)
        self.assertEqual([(r.start, r.end) for r in regions],
                         [(401, 800), (1401, 1800)])

    def test_sequence_exactly_one_window_long(self):
        regions = self.detect(island(200))
        self.assertEqual([(r.start, r.end) for r in regions], [(1, 200)])

    def test_island_at_very_end_found_via_end_flush_window(self):
        # n = 400, window 200, step 300: grid gives only [0,200) (flank).
        # Without the end-flush window the island at 201-400 would be missed.
        regions = self.detect(flank(200) + island(200), step=300)
        self.assertEqual([(r.start, r.end) for r in regions], [(201, 400)])

    def test_island_at_very_start(self):
        regions = self.detect(island(400) + flank(400))
        self.assertEqual([(r.start, r.end) for r in regions], [(1, 400)])

    def test_window_larger_than_sequence_scans_whole_sequence(self):
        regions = self.detect(island(248), window=300)
        self.assertEqual([(r.start, r.end) for r in regions], [(1, 248)])

    def test_sequence_shorter_than_min_length_finds_nothing(self):
        self.assertEqual(self.detect(island(192)), [])

    def test_empty_sequence_finds_nothing(self):
        self.assertEqual(self.detect(""), [])

    def test_gc_rich_low_obs_exp_sequence_finds_nothing(self):
        self.assertEqual(self.detect("G" * 300 + "C" * 300), [])

    def test_at_rich_sequence_finds_nothing(self):
        self.assertEqual(self.detect(flank(1000)), [])

    def test_lowercase_input_gives_same_regions(self):
        sequence = flank(400) + island(400) + flank(400)
        self.assertEqual(self.detect(sequence.lower()), self.detect(sequence))

    def test_n_gap_inside_island_splits_region_because_windows_ineligible(self):
        # island 401-800 with N at 561-600 (40 N). Any 200-bp window that
        # includes them has 20 % N > 5 %, so those windows are ineligible;
        # step 100 windows [400,600) and [500,700) both include the N run.
        # Only [600,800) survives -> 601-800.
        sequence = (flank(400) + island(160) + "N" * 40 + island(200)
                    + flank(400))
        regions = self.detect(sequence)
        self.assertEqual([(r.start, r.end) for r in regions], [(601, 800)])

    def test_n_gap_is_tolerated_when_fraction_allowed(self):
        sequence = (flank(400) + island(160) + "N" * 40 + island(200)
                    + flank(400))
        regions = self.detect(sequence, max_ambiguous_fraction=0.5)
        self.assertEqual(len(regions), 1)
        self.assertEqual(regions[0].metrics.ambiguous_count, 40)

    def test_profile_reports_every_window(self):
        sequence = flank(400) + island(400) + flank(400)
        params = cpg.DetectionParameters(window_size_bp=200, step_size_bp=100)
        record = cpg.SequenceRecord("t", "", sequence)
        result = cpg.analyse_record(record, params, include_profile=True)
        starts = [point.start0 for point in result.profile]
        self.assertEqual(starts, list(range(0, 1001, 100)))
        # Window [400,600) lies inside the island: GC = 75 %.
        self.assertAlmostEqual(result.profile[4].gc_percent, 75.0)
        # Window [0,200) is pure AT: GC = 0 %, Obs/Exp undefined (C = G = 0).
        self.assertEqual(result.profile[0].gc_percent, 0.0)
        self.assertIsNone(result.profile[0].obs_exp_cpg)

    def test_coordinate_conversion(self):
        self.assertEqual(cpg.to_one_based_inclusive(0, 10), (1, 10))
        self.assertEqual(cpg.to_one_based_inclusive(400, 800), (401, 800))
        sequence = flank(400) + island(400) + flank(400)
        region = self.detect(sequence)[0]
        # Round trip back to a Python slice must give exactly the island.
        self.assertEqual(sequence[region.start - 1:region.end], island(400))


class TestSyntheticValidation(unittest.TestCase):
    """Validate against data/sample.fasta (synthetic, known structure)."""

    TRUTH = {
        "synthetic_island_01": [(401, 900)],
        "synthetic_two_islands_02": [(301, 650), (1151, 1750)],
        "synthetic_gc_rich_cpg_depleted_03": [],
        "synthetic_at_rich_background_04": [],
        "synthetic_short_cpg_rich_05": [],
        "synthetic_island_with_Ns_06": [(301, 700)],
    }

    @classmethod
    def setUpClass(cls):
        records = cpg.read_fasta_file(str(PROJECT_ROOT / "data" /
                                          "sample.fasta"))
        cls.results = {r.identifier: cpg.analyse_record(r, DEFAULTS)
                       for r in records}

    def test_all_expected_records_present(self):
        self.assertEqual(set(self.results), set(self.TRUTH))

    def test_each_true_island_is_recovered_with_bounded_boundary_error(self):
        # Window-based boundaries are approximate: allow one window length
        # of error. Every true island must be overlapped by a predicted
        # region covering >= 80 % of it, and no predicted region may stray
        # more than one window length from a true island.
        slack = DEFAULTS.window_size_bp
        for identifier, truth in self.TRUTH.items():
            if not truth:
                continue
            predicted = [(r.start, r.end)
                         for r in self.results[identifier].regions]
            for true_start, true_end in truth:
                # Predicted regions never overlap each other, so coverage of
                # the true island is the sum of the individual overlaps.
                covered = sum(max(0, min(end, true_end)
                                  - max(start, true_start) + 1)
                              for start, end in predicted)
                self.assertGreaterEqual(
                    covered / (true_end - true_start + 1), 0.8,
                    msg=f"{identifier}: island {true_start}-{true_end} "
                        f"poorly covered by {predicted}")
            for start, end in predicted:
                self.assertTrue(
                    any(start >= ts - slack and end <= te + slack
                        for ts, te in truth),
                    msg=f"{identifier}: region {start}-{end} too far from "
                        f"truth {truth}")

    def test_negative_controls_yield_no_regions(self):
        for identifier, truth in self.TRUTH.items():
            if not truth:
                self.assertEqual(self.results[identifier].regions, [],
                                 msg=identifier)

    def test_number_of_regions_matches_truth_except_documented_n_gap(self):
        for identifier, truth in self.TRUTH.items():
            if identifier == "synthetic_island_with_Ns_06":
                continue  # the 30-N gap splits the island; tested below
            self.assertEqual(len(self.results[identifier].regions),
                             len(truth), msg=identifier)

    def test_n_gap_record_is_split_at_the_gap(self):
        regions = self.results["synthetic_island_with_Ns_06"].regions
        self.assertEqual(len(regions), 2)
        # The split happens at the N run (451-480): the first piece ends
        # inside/before it and the second starts after the run's start.
        self.assertLess(regions[0].end, 451 + 30)
        self.assertGreater(regions[1].start, 450)


class TestExport(unittest.TestCase):
    def setUp(self):
        params = cpg.DetectionParameters(window_size_bp=200, step_size_bp=100)
        sequence = flank(400) + island(400) + flank(400)
        self.params = params
        self.results = [
            cpg.analyse_record(cpg.SequenceRecord("one", "", sequence),
                               params),
            cpg.analyse_record(cpg.SequenceRecord("two", "", flank(300)),
                               params),
        ]

    def test_csv_export_content(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "sub" / "out.csv"
            cpg.export_results(self.results, self.params, str(path))
            with open(path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual((row["sequence_id"], row["start_bp"], row["end_bp"],
                          row["length_bp"]), ("one", "401", "800", "400"))
        self.assertEqual(row["gc_percent"], "75.00")
        self.assertEqual(row["obs_exp_cpg"], "2.667")
        self.assertEqual(row["cpg_count"], "150")

    def test_json_export_content(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "out.json"
            cpg.export_results(self.results, self.params, str(path))
            document = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(document["coordinate_system"], "1-based, inclusive")
        self.assertEqual(len(document["records"]), 2)
        self.assertEqual(document["records"][0]["regions"][0]["start_bp"], 401)
        self.assertEqual(document["records"][1]["regions"], [])

    def test_svg_is_well_formed_xml(self):
        import xml.etree.ElementTree as element_tree
        params = self.params
        record = cpg.SequenceRecord("a<b>&c", "",
                                    flank(400) + island(400) + flank(400))
        result = cpg.analyse_record(record, params, include_profile=True)
        root = element_tree.fromstring(cpg.build_svg_profile(result, params))
        self.assertTrue(root.tag.endswith("svg"))

    def test_text_map_marks_expected_columns(self):
        # 100 bp, region 26-75 -> with width 4 (25 bp per cell): .##.
        region = cpg.CandidateRegion("x", 26, 75,
                                     cpg.compute_metrics("CG" * 25))
        self.assertEqual(cpg.render_text_map(100, [region], width=4), ".##.")

    def test_safe_file_stem(self):
        self.assertEqual(cpg.safe_file_stem("chr1:100-200 |x"),
                         "chr1_100-200_x")
        self.assertEqual(cpg.safe_file_stem("///"), "record")


class TestCommandLine(unittest.TestCase):
    @staticmethod
    def run_cli(*argv):
        stdout, stderr = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(stdout), \
                contextlib.redirect_stderr(stderr):
            try:
                code = cpg.main(list(argv))
            except SystemExit as exit_request:
                code = exit_request.code
        return code, stdout.getvalue(), stderr.getvalue()

    def test_direct_sequence_success(self):
        code, out, _ = self.run_cli("-s", island(400), "--window-size", "200",
                                    "--step-size", "100")
        self.assertEqual(code, 0)
        self.assertIn("Candidate regions: 1", out)
        table_rows = [line.split() for line in out.splitlines()
                      if line.strip().startswith("1 ")]
        # Row layout: number, start, end, length, GC%, Obs/Exp, CpG.
        self.assertEqual(table_rows[0][:4], ["1", "1", "400", "400"])

    def test_missing_file_gives_clear_error(self):
        code, _, err = self.run_cli("-f", "no_such_file.fasta")
        self.assertEqual(code, 1)
        self.assertIn("file not found: no_such_file.fasta", err)

    def test_invalid_sequence_gives_clear_error(self):
        code, _, err = self.run_cli("-s", "ACGTX")
        self.assertEqual(code, 1)
        self.assertIn("'X' at position 5", err)

    def test_ambiguous_reject_flag(self):
        code, _, err = self.run_cli("-s", "ACGTN", "--ambiguous", "reject")
        self.assertEqual(code, 1)
        self.assertIn("ambiguous", err)

    def test_invalid_parameter_exits_with_usage_error(self):
        code, _, err = self.run_cli("-s", "ACGT", "--step-size", "0")
        self.assertEqual(code, 2)
        self.assertIn("step-size must be at least 1", err)

    def test_empty_sequence_is_reported_not_crashed(self):
        code, out, _ = self.run_cli("-s", "")
        self.assertEqual(code, 0)
        self.assertIn("Candidate regions: 0", out)
        self.assertIn("GC NA %", out)

    def test_no_input_source_is_a_usage_error(self):
        code, _, _ = self.run_cli()
        self.assertEqual(code, 2)

    def test_output_and_plot_files_are_written(self):
        sample = str(PROJECT_ROOT / "data" / "sample.fasta")
        with tempfile.TemporaryDirectory() as folder:
            csv_path = str(Path(folder) / "regions.csv")
            plot_dir = str(Path(folder) / "plots")
            code, out, _ = self.run_cli("-f", sample, "-o", csv_path,
                                        "--plot-dir", plot_dir)
            self.assertEqual(code, 0)
            self.assertTrue(Path(csv_path).is_file())
            self.assertEqual(len(list(Path(plot_dir).glob("*.svg"))), 6)

    def test_output_is_deterministic(self):
        sample = str(PROJECT_ROOT / "data" / "sample.fasta")
        first = self.run_cli("-f", sample)
        second = self.run_cli("-f", sample)
        self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main(verbosity=2)
