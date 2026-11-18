#!/usr/bin/env python3
"""Generate data/sample.fasta: SYNTHETIC sequences with known structure.

These sequences are simulated, not real genomic DNA. They mimic the broad
compositional contrast between vertebrate "bulk" DNA (AT-rich, CpG-depleted)
and CpG-island-like DNA (GC-rich, little CpG depletion), so the detector can
be validated against a known answer. A tiny linear-congruential generator is
used instead of the ``random`` module so the output is identical on every
Python version and platform.

Usage:  python scripts/generate_sample_fasta.py
"""

from pathlib import Path
from typing import Dict, List

LCG_MULTIPLIER = 6364136223846793005
LCG_INCREMENT = 1442695040888963407
LCG_MODULUS = 2 ** 64
FASTA_LINE_WIDTH = 70
SEED = 20240611

# Base probabilities and the factor applied to P(G) right after a C.
# A factor of 1.0 means no CpG depletion; small factors mean strong depletion.
ISLAND_LIKE = {"A": 0.20, "C": 0.30, "G": 0.30, "T": 0.20, "g_after_c": 0.80}
BULK_LIKE = {"A": 0.31, "C": 0.19, "G": 0.19, "T": 0.31, "g_after_c": 0.20}
GC_RICH_DEPLETED = {"A": 0.17, "C": 0.33, "G": 0.33, "T": 0.17,
                    "g_after_c": 0.15}


class DeterministicRandom:
    """64-bit LCG producing floats in [0, 1)."""

    def __init__(self, seed: int) -> None:
        self.state = seed % LCG_MODULUS

    def next_float(self) -> float:
        self.state = (self.state * LCG_MULTIPLIER + LCG_INCREMENT) % LCG_MODULUS
        return (self.state >> 11) / float(2 ** 53)


def generate_segment(rng: DeterministicRandom, length: int,
                     model: Dict[str, float]) -> str:
    """Generate a segment base by base from the given composition model."""
    bases: List[str] = []
    for _ in range(length):
        weights = {base: model[base] for base in "ACGT"}
        if bases and bases[-1] == "C":
            weights["G"] *= model["g_after_c"]
        total = sum(weights.values())
        draw = rng.next_float() * total
        cumulative = 0.0
        chosen = "T"
        for base in "ACGT":
            cumulative += weights[base]
            if draw < cumulative:
                chosen = base
                break
        bases.append(chosen)
    return "".join(bases)


def wrap(sequence: str) -> str:
    return "\n".join(sequence[i:i + FASTA_LINE_WIDTH]
                     for i in range(0, len(sequence), FASTA_LINE_WIDTH))


def main() -> None:
    rng = DeterministicRandom(SEED)
    records = []

    # 1. One island in the middle: truth 401-900.
    sequence = (generate_segment(rng, 400, BULK_LIKE)
                + generate_segment(rng, 500, ISLAND_LIKE)
                + generate_segment(rng, 400, BULK_LIKE))
    records.append(("synthetic_island_01",
                    "SYNTHETIC; 1300 bp; truth: island at 401-900", sequence))

    # 2. Two islands: truth 301-650 and 1151-1750.
    sequence = (generate_segment(rng, 300, BULK_LIKE)
                + generate_segment(rng, 350, ISLAND_LIKE)
                + generate_segment(rng, 500, BULK_LIKE)
                + generate_segment(rng, 600, ISLAND_LIKE)
                + generate_segment(rng, 300, BULK_LIKE))
    records.append(("synthetic_two_islands_02",
                    "SYNTHETIC; 2050 bp; truth: islands at 301-650 and "
                    "1151-1750", sequence))

    # 3. GC-rich decoy: high GC but strongly CpG-depleted -> no island.
    records.append(("synthetic_gc_rich_cpg_depleted_03",
                    "SYNTHETIC; 1000 bp; truth: no island (Obs/Exp low)",
                    generate_segment(rng, 1000, GC_RICH_DEPLETED)))

    # 4. AT-rich bulk background -> no island.
    records.append(("synthetic_at_rich_background_04",
                    "SYNTHETIC; 800 bp; truth: no island",
                    generate_segment(rng, 800, BULK_LIKE)))

    # 5. CpG-rich but only 150 bp -> too short.
    records.append(("synthetic_short_cpg_rich_05",
                    "SYNTHETIC; 150 bp; truth: no island (too short)",
                    generate_segment(rng, 150, ISLAND_LIKE)))

    # 6. Island containing a 30-base run of N (assembly gap): truth 301-700,
    #    N run at 451-480.
    island = generate_segment(rng, 400, ISLAND_LIKE)
    island = island[:150] + "N" * 30 + island[180:]
    sequence = (generate_segment(rng, 300, BULK_LIKE) + island
                + generate_segment(rng, 300, BULK_LIKE))
    records.append(("synthetic_island_with_Ns_06",
                    "SYNTHETIC; 1000 bp; truth: island at 301-700 with "
                    "N gap at 451-480", sequence))

    output_path = Path(__file__).resolve().parent.parent / "data" / \
        "sample.fasta"
    with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
        for identifier, description, sequence in records:
            handle.write(f">{identifier} {description}\n{wrap(sequence)}\n")
    print(f"Wrote {len(records)} records to {output_path}")


if __name__ == "__main__":
    main()
