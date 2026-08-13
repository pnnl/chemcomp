#!/usr/bin/env python3
"""Convert MotifsDB CSV rows into TableGen motif records.

Used to (re)generate curated fixture records for
include/ChemComp/Conversion/CRNToDB/MotifsFixture.td. Run offline; its output
is meant to be reviewed and pasted, not wired into the build.

Emitting all 43k motifs as TableGen is deliberately not supported: parsing that
many records would dominate build time, and the CSVs are untracked so the build
must not depend on them. Use csv-to-json.py for the full database and keep the
compiled-in fixture small.

    ./chemcomp/utils/csv-to-td.py \
        --reactions MotifsDB/reactions_reduced.csv \
        --compounds MotifsDB/compounds_reduced.csv \
        --id rxn00001 --id rxn00006 \
        --slice-name PaperTable1
"""

import argparse
import re
import sys

from motifsdb_parse import parse_database


def record_name(motif_id: str) -> str:
    """TableGen def name for a motif id: rxn00001 -> Rxn00001."""
    cleaned = re.sub(r"[^A-Za-z0-9_]", "_", motif_id)
    return cleaned[:1].upper() + cleaned[1:]


def tblgen_string(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def emit_motif(motif, out) -> None:
    species = ", ".join(tblgen_string(s) for s in motif.species)
    stoich = ", ".join(str(v) for v in motif.stoich)

    print(f"def {record_name(motif.motif_id)} : Motif<"
          f"{tblgen_string(motif.motif_id)},", file=out)
    print(f"    [{species}], [{stoich}]> {{", file=out)
    # Real-valued fields are strings because TableGen has no float type; the
    # empty default already means unknown, so omit rather than write "".
    if motif.delta_g is not None:
        print(f'  let DeltaG = "{motif.delta_g}";', file=out)
    if motif.delta_g_err is not None:
        print(f'  let DeltaGErr = "{motif.delta_g_err}";', file=out)
    print(f"  let Source = {tblgen_string(motif.source)};", file=out)
    print("}", file=out)
    print(file=out)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reactions", required=True)
    ap.add_argument("--compounds", required=True)
    ap.add_argument(
        "--id",
        action="append",
        default=[],
        dest="ids",
        required=True,
        help="motif id to emit (repeatable)",
    )
    ap.add_argument(
        "--slice-name",
        default="",
        help="also emit a Slice def containing the selected motifs",
    )
    ap.add_argument("-o", "--output", default="-")
    args = ap.parse_args()

    motifs, _ = parse_database(args.reactions, args.compounds)
    by_id = {m.motif_id: m for m in motifs}

    missing = [i for i in args.ids if i not in by_id]
    if missing:
        # A requested id may have been dropped by the parser (non-integer
        # coefficients, one-sided reaction); say so rather than emit silence.
        print(f"error: not found or skipped during parsing: {', '.join(missing)}",
              file=sys.stderr)
        return 1

    selected = [by_id[i] for i in sorted(args.ids)]
    out = sys.stdout if args.output == "-" else open(args.output, "w",
                                                     encoding="utf-8")
    try:
        for motif in selected:
            emit_motif(motif, out)

        if args.slice_name:
            members = ", ".join(record_name(m.motif_id) for m in selected)
            print(f"def {args.slice_name} : Slice<"
                  f'"{args.slice_name}", [{members}]>;', file=out)
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
