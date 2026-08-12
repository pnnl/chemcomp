"""Parse the ModelSEED-derived MotifsDB CSVs into motif records.

Shared by csv-to-td.py and csv-to-json.py so there is exactly one parser and
the TableGen fixture cannot drift from the runtime JSON database.

These scripts are run offline by a developer, never as a build step: the CSVs
are untracked in the MotifsDB submodule, so a fresh clone does not have them
and a build-time dependency would break for everyone else.

Data hazards handled here, all measured against the real files
(43,774 reactions / 33,992 compounds):

  * 5,272 compound names are quoted and contain commas, so a real CSV reader
    is mandatory -- splitting on ',' corrupts them.
  * `10000000` is a sentinel meaning "no thermodynamic data", not a value. It
    appears in 18,851 of 43,774 rows (43.1%) and becomes None here.
  * 1,397 coefficients are non-integer (decimals, and scientific notation such
    as `6.06e-05`), affecting ~295 reactions. The stoich schema is integral, so
    those reactions are skipped and reported.
  * Arrows are `<=>` (27,676), `=>` (12,767) and `<=` (3,331). `<=>` contains
    both others as substrings, so it must be tested first; `<=` is a reverse
    reaction whose sides are swapped to normalize to reactant-first.
  * 140 compound names are shared by different cpd ids, and 2 are empty. Names
    are disambiguated by appending the id when ambiguous.
  * Compartment tags are `[0]` (194,165), `[1]` (7,136) and `[2]` (97). The
    schema has no compartment field, so non-zero compartments are folded into
    the species name to keep them from colliding with the cytosolic form.
  * 11 reactions are empty and 76 have an empty side; both are skipped.
  * Neither file ends with a trailing newline.
"""

import csv
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# Sentinel used by ModelSEED for missing thermodynamic data.
UNKNOWN_DELTAG = 1e7

# `(coeff) cpdNNNNN[compartment]`
_TERM_RE = re.compile(r"\(([^)]*)\)\s+(cpd\d+)\[(\d+)\]")

# Longest first: `<=>` contains both `=>` and `<=`.
_ARROWS = ("<=>", "=>", "<=")


@dataclass
class Motif:
    motif_id: str
    species: List[str]
    stoich: List[int]
    delta_g: Optional[float]
    delta_g_err: Optional[float]
    source: str


@dataclass
class ParseStats:
    total_rows: int = 0
    emitted: int = 0
    skipped_no_arrow: int = 0
    skipped_empty_side: int = 0
    skipped_non_integer: int = 0
    skipped_duplicate_species: int = 0
    reasons: Dict[str, int] = field(default_factory=dict)

    def note(self, reason: str) -> None:
        self.reasons[reason] = self.reasons.get(reason, 0) + 1


def load_compounds(path: str) -> Dict[str, str]:
    """Map cpd id -> display name, disambiguating names shared by several ids.

    The target schema keys on names rather than ids, so a name used by more
    than one compound would silently merge two distinct chemicals. Those get
    the id appended; unique names are left alone so the common case stays
    readable.
    """
    raw: Dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            cid = row["id"].strip()
            name = (row.get("name") or "").strip()
            raw[cid] = name or cid

    counts: Dict[str, int] = {}
    for name in raw.values():
        counts[name] = counts.get(name, 0) + 1

    return {
        cid: (f"{name} ({cid})" if counts.get(name, 0) > 1 else name)
        for cid, name in raw.items()
    }


def _species_key(names: Dict[str, str], cpd: str, compartment: str) -> str:
    """Name for a compound in a compartment.

    Compartment 0 (cytosol) is the overwhelming majority and stays bare; other
    compartments are suffixed so the same compound in two compartments does not
    collapse into one species.
    """
    base = names.get(cpd, cpd)
    return base if compartment == "0" else f"{base}[{compartment}]"


def _parse_optional_float(text: str) -> Optional[float]:
    try:
        value = float(text)
    except (TypeError, ValueError):
        return None
    return None if value >= UNKNOWN_DELTAG else value


def _split_equation(equation: str) -> Optional[Tuple[str, str]]:
    """Split into (reactants, products), normalizing reverse reactions."""
    for arrow in _ARROWS:
        if arrow in equation:
            left, right = equation.split(arrow, 1)
            # `<=` runs right-to-left, so swap to keep reactants first.
            return (right, left) if arrow == "<=" else (left, right)
    return None


def parse_reaction(
    row: Dict[str, str], names: Dict[str, str], stats: ParseStats
) -> Optional[Motif]:
    motif_id = row["motif_id"].strip()

    sides = _split_equation(row["equation"])
    if sides is None:
        stats.skipped_no_arrow += 1
        stats.note("no arrow")
        return None

    lhs, rhs = sides
    species: List[str] = []
    stoich: List[int] = []

    for side, sign in ((lhs, -1), (rhs, 1)):
        for coeff_text, cpd, compartment in _TERM_RE.findall(side):
            try:
                value = float(coeff_text)
            except ValueError:
                stats.skipped_non_integer += 1
                stats.note(f"unparseable coefficient in {motif_id}")
                return None
            if value != int(value):
                # Biomass-style reactions carry fractional coefficients; the
                # stoich schema and stoich_matrix are integral.
                stats.skipped_non_integer += 1
                stats.note(f"non-integer coefficient in {motif_id}")
                return None
            species.append(_species_key(names, cpd, compartment))
            stoich.append(sign * int(value))

    if not species:
        stats.skipped_empty_side += 1
        stats.note("empty reaction")
        return None

    # A species on both sides appears twice; combine into one net coefficient
    # rather than rejecting, since that is a real (and correct) reading.
    combined: Dict[str, int] = {}
    order: List[str] = []
    for name, coeff in zip(species, stoich):
        if name not in combined:
            order.append(name)
            combined[name] = 0
        combined[name] += coeff

    species = [n for n in order if combined[n] != 0]
    stoich = [combined[n] for n in species]

    if not species:
        # Everything cancelled: not a reaction.
        stats.skipped_empty_side += 1
        stats.note("fully cancelling reaction")
        return None

    if not any(v < 0 for v in stoich) or not any(v > 0 for v in stoich):
        # Exchange/sink reactions have only one side; they cannot be matched
        # against a CRN reaction that consumes and produces.
        stats.skipped_empty_side += 1
        stats.note("one-sided reaction")
        return None

    return Motif(
        motif_id=motif_id,
        species=species,
        stoich=stoich,
        delta_g=_parse_optional_float(row.get("deltag", "")),
        delta_g_err=_parse_optional_float(row.get("deltagerr", "")),
        # The CSVs carry no cross-reference column -- `grep -c KEGG` is 0 on
        # both -- so only the ModelSEED accession is derivable locally. The
        # paper's `KEGG:...|ModelSEED:...` strings were hand-authored for its
        # ten-record sample.
        source=f"ModelSEED:{motif_id}",
    )


def parse_database(
    reactions_path: str, compounds_path: str
) -> Tuple[List[Motif], ParseStats]:
    names = load_compounds(compounds_path)
    stats = ParseStats()
    motifs: List[Motif] = []

    with open(reactions_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            stats.total_rows += 1
            motif = parse_reaction(row, names, stats)
            if motif is not None:
                motifs.append(motif)
                stats.emitted += 1

    motifs.sort(key=lambda m: m.motif_id)
    return motifs, stats


def format_stats(stats: ParseStats) -> str:
    lines = [
        f"rows read:            {stats.total_rows}",
        f"motifs emitted:       {stats.emitted}",
        f"skipped (no arrow):   {stats.skipped_no_arrow}",
        f"skipped (empty side): {stats.skipped_empty_side}",
        f"skipped (non-integer):{stats.skipped_non_integer}",
    ]
    return "\n".join(lines)
