#!/usr/bin/env python3
"""Convert the MotifsDB CSVs into the runtime JSON motif database.

This is the full-scale database the crn-db-match pass loads via `db=<path>`.
Run offline; the output is not committed and not a build dependency.

    ./chemcomp/utils/csv-to-json.py \
        --reactions MotifsDB/reactions_reduced.csv \
        --compounds MotifsDB/compounds_reduced.csv \
        -o /tmp/motifs.json

Slicing options exist because the full database does not bound the matcher's
candidate sets: the single signature (-1,-1,1,1) holds 11,421 of 43,455 motifs,
so a species-restricted or capped subset is usually what you want.
"""

import argparse
import json
import sys

from motifsdb_parse import format_stats, parse_database


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--compounds", required=True)
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument(
        "--require-species",
        action="append",
        default=[],
        help="keep only motifs containing this species (repeatable)",
    )
    ap.add_argument(
        "--max-motifs",
        type=int,
        default=0,
        help="keep at most this many motifs (0 = unlimited)",
    )
    ap.add_argument(
        "--indent",
        type=int,
        default=1,
        help="JSON indentation; 0 emits one compact line per file",
    )
    args = ap.parse_args()

    motifs, stats = parse_database(args.reactions, args.compounds)

    if args.require_species:
        wanted = set(args.require_species)
        motifs = [m for m in motifs if wanted & set(m.species)]

    if args.max_motifs:
        motifs = motifs[: args.max_motifs]

    records = []
    for m in motifs:
        record = {
            "motif_id": m.motif_id,
            "species": m.species,
            "stoich": m.stoich,
            "source": m.source,
        }
        # Absent rather than null when unknown: the loader treats a missing key
        # and an explicit null identically, and omitting keeps the file smaller.
        if m.delta_g is not None:
            record["deltaG"] = m.delta_g
        if m.delta_g_err is not None:
            record["deltaGerr"] = m.delta_g_err
        records.append(record)

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(records, handle, indent=args.indent or None)
        handle.write("\n")

    print(format_stats(stats), file=sys.stderr)
    print(f"written:              {len(records)} -> {args.output}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
