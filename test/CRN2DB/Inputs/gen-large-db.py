#!/usr/bin/env python3
"""Synthesize a motif database with the real database's signature skew.

The scale test needs a database large enough that candidate truncation and the
state-count guard actually engage, but MotifsDB's CSVs are untracked in the
submodule, so a fresh clone cannot generate one from real data. This produces a
stand-in with the property that matters -- a few signatures holding most of the
motifs -- without any external input.

The skew reproduced here was measured over the 43,455 integer-coefficient
motifs of the real database:

    (-1,-1,1,1)     11,421   26.3%
    (-1,-1,1,1,1)    6,999   16.1%
    (-1,1)           5,153   11.9%
    (-1,-1,-1,1,1)   2,554    5.9%
    everything else         39.8%  spread over ~909 signatures

Regenerate with:

    ./gen-large-db.py -o large-db.json --motifs 5000
"""

import argparse
import json

# (signature, share of the database) for the head of the real distribution.
SKEW = [
    ((-1, -1, 1, 1), 0.263),
    ((-1, -1, 1, 1, 1), 0.161),
    ((-1, 1), 0.119),
    ((-1, -1, -1, 1, 1), 0.059),
]

# The tail: distinct signatures splitting the remaining share evenly.
TAIL = [
    (-1, -1, -1, 1, 1, 1, 1),
    (-1, 1, 1),
    (-1, -1, 1),
    (-2, -1, 1, 2),
    (-1, 1, 2),
    (-2, 1, 2),
    (-1, -1, 2),
    (-3, 1, 1, 1),
]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--motifs", type=int, default=5000)
    args = ap.parse_args()

    records = []

    def emit(signature):
        index = len(records)
        # Species names are drawn from a pool much smaller than the motif count
        # so motifs genuinely share species, which is what makes AC-3 and the
        # injectivity constraints do real work.
        species = [f"C{(index * 7 + i * 13) % 512:04d}" for i in
                   range(len(signature))]
        # A motif repeating a species is rejected by the loader; the stride
        # above can collide, so fall back to unique names when it does.
        if len(set(species)) != len(species):
            species = [f"C{index:05d}_{i}" for i in range(len(signature))]
        records.append({
            "motif_id": f"syn{index:06d}",
            "species": species,
            "stoich": list(signature),
            "source": "Synthetic:scale-test",
        })

    for signature, share in SKEW:
        for _ in range(int(args.motifs * share)):
            emit(signature)

    while len(records) < args.motifs:
        emit(TAIL[len(records) % len(TAIL)])

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(records, handle)
        handle.write("\n")

    print(f"wrote {len(records)} motifs to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
