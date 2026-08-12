# CRN2DB: matching abstract CRNs to a reaction database

ChemComp compiles ODEs into *abstract* CRNs — `acr.crn_generic` operations whose
species are anonymous column indices of a `stoich_matrix`. Running that CRN on a
real biochemical device needs each abstract species bound to a concrete chemical
and each abstract reaction bound to a real reaction with compatible
stoichiometry. The `crn-db-match` pass does that binding.

Implements *"CRN2DB: Chemical Reaction Networks to Database Matching via Dynamic
Programming in ChemComp"* (`papers/cf26-crn2db`). Where this implementation
departs from the paper, the reasons are recorded at the end — several are
corrections to errors in the paper itself.

## Using it

```bash
# Compiled-in fixture, no external data needed.
standalone-opt -crn-db-match='db-slice=OscillatorSlice' input.mlir

# Full database from JSON.
standalone-opt -crn-db-match='db=motifs.json' input.mlir

# The whole flow, out to the emulator.
standalone-opt -crn-db-match='db-slice=OscillatorSlice' input.mlir \
  | standalone-translate -mlir-to-json
```

The pass reads `stoich_matrix` and annotates the same operation in place with
two discardable attributes:

```mlir
%out = acr.crn_generic {
    stoich_matrix = dense<...> : tensor<4x6xi32>,
    species_names = ["ATP", "ADP", "AMP", "Phosphate", "H2O", "H+"],
    motif_ids = ["rxnosc01", "rxnosc02", "rxnosc03", "rxnosc04"]
  } (%species, %rates) : ...
```

Unmatched reactions and unbound species get the empty string, so both arrays
stay aligned with the matrix. Partial coverage is a warning, not an error: the
paper's flow expects downstream passes to consume partially instantiated CRNs.

## Options

| Option | Default | Meaning |
|---|---|---|
| `db` | `""` | JSON database path; empty selects the compiled-in fixture |
| `db-slice` | `""` | Named fixture slice; empty means every fixture motif |
| `filter-species` | — | Keep only motifs containing at least one of these |
| `filter-source` | `""` | Keep only motifs whose source has this prefix |
| `max-candidates` | 256 | Cap on candidates per reaction; **truncation emits a remark** |
| `max-per-signature` | 0 | Cap on motifs loaded per signature; shrinks the database, silent |
| `wc` / `ws` / `wr` | 2.0 / 1.0 / 0.5 | Coverage, stoichiometry and rate weights |
| `sigma` | 1.0 | Width of the log-space Gaussian in `S_r` |
| `beam-width` | 0 | 0 runs the exact DP; positive prunes to that many states |
| `strict-signature` | false | Exact signature equality instead of superset lookup |
| `global-injectivity` | false | Enforce a globally injective assignment |
| `dump-db` | false | Print the loaded database and stop |
| `dump-preprocessing` | false | Print ordering, frontiers, domains and candidates |

`max-candidates` and `max-per-signature` are different mechanisms and it matters
which you reach for. The first caps what the matcher considers and reports what
it dropped; the second discards motifs at load time, so they are never
candidates and nothing is reported. Capping at load time would silently suppress
the remark that tells you the result stopped being exact.

## Database

Two sources, one in-memory representation, so the matcher cannot tell them
apart. `test/CRN2DB/db-loader-equivalence.mlir` asserts that by running one set
of expectations against both.

**Compiled-in fixture** — ~50 curated motifs in
`include/Standalone/CRN2DB/MotifsFixture.td`, turned into a static table by
`chemcomp-tblgen -gen-motif-db`. Tests use this so `check-standalone` needs no
external file.

**Runtime JSON** — the full database, generated offline:

```bash
./utils/csv-to-json.py --reactions MotifsDB/reactions_reduced.csv \
                       --compounds MotifsDB/compounds_reduced.csv \
                       -o motifs.json
```

Only the fixture is compiled in. The full database is 43,390 motifs; parsing
that as TableGen would dominate build time, and the MotifsDB CSVs are untracked
in the submodule, so a build-time dependency would break any fresh clone.

### Extending the schema

Typed properties live on `Motif` subclasses, so a new column costs a subclass
and one line in `chemcomp-tblgen/MotifDBGen.cpp` — nothing downstream changes.
`SolubleMotif` is the worked example:

```tablegen
class SolubleMotif<string id, list<string> sp, list<int> st> : Motif<id, sp, st> {
  bit Soluble = 1;
}

def SolubleSlice : Slice<"SolubleSlice", [...]>;
```

TableGen has no float type, so real-valued fields are decimal strings parsed by
the backend, and no set or map type, so species and coefficients are parallel
lists that the generated C++ turns into a sorted set and a map. The empty string
means unknown, which stays distinct from any numeric value.

## Algorithm

1. **Species-sharing graph** — species are nodes, two adjacent when a reaction
   touches both.
2. **Min-fill elimination ordering** over species, then a derived reaction order.
3. **Candidate collection** by stoichiometric signature — the sorted multiset of
   non-zero coefficients.
4. **AC-3 arc consistency** pruning each species' domain of concrete chemicals.
5. **Frontier DP** over the reaction order, keeping only species that are both
   already seen and still to come.

Inspect stages 1–4 with `dump-preprocessing`.

## Deviations from the paper

### Errors found in the paper

**`|Cand| <= 100` is false.** Section 4.3 claims signature indexing bounds the
candidate set. Measured over the real database: 43,455 integer-coefficient
motifs across only 913 distinct signatures, with `(-1,-1,1,1)` alone holding
11,421 — and 90.2% of motifs in buckets larger than 100. The coverage-correct
superset lookup returns **30,649 candidates**, 70% of the database, for that one
shape. The stated complexity bound does not hold at full scale; `max-candidates`
is what makes matching tractable, and it makes the result heuristic, hence the
remark.

**The oscillator's treewidth is 3, not 2.** Section 5.4 claims `W = 2` for its
example, reasoning that "species A, B, D, E each appear in two reactions". But
section 4.1 defines frontier width as species simultaneously *live*, not
occurrences per species. Exhaustive search over all 24 reaction orderings gives a
minimum achievable maximum frontier width of 3. Our min-fill ordering achieves
that true optimum. See `test/CRN2DB/preprocessing-oscillator.mlir`.

**Table 1's `rxn00003` is transcribed without flipping its arrow.** The source
equation uses `<=`, a reverse reaction, so its reactants are Pyruvate and H+ —
not the products the paper lists. The fixture uses the CSV-derived form.

### Gaps the paper leaves open

**Species order to reaction order.** Min-fill runs over species but the DP
consumes a reaction ordering, and the paper does not say how one becomes the
other. A reaction is fully determined only once every species it touches has been
introduced, so reactions are keyed by the latest elimination position in their
support, with deterministic tie-breaks.

**Global injectivity is not representable in the paper's DP state.** Section 3.3
requires an injective `phi`, but `trim()` drops a binding once its species leaves
the frontier, after which the chemical can be rebound. The default is
paper-faithful frontier-local injectivity **plus a warning** naming the reused
chemicals. `global-injectivity=1` carries retired chemicals in the state key and
is genuinely injective, at the cost of a state space exponential in the species
count rather than the treewidth. `test/CRN2DB/dp-injectivity.mlir` runs one input
through both modes.

**Exact-signature indexing contradicts the coverage rule.** Section 3.2 permits
`phi(supp(r)) ⊆ D_j`, so a motif with extra species legitimately matches — but it
has a longer signature and exact-key lookup misses it. The default is a superset
lookup; `strict-signature=1` restores paper-literal behaviour.

**No tie-breaking rule is given**, yet results are compared by FileCheck. Imposed
order: score, then match count, then motif-id sequence, then binding. Every
observable iteration runs over sorted vectors, never a raw hash map.

**`kappa_j` has no scalar field** in the schema, only `rate_range`. The nominal
rate is the geometric mean of that range, since `S_r` compares rates in log
space. Absent range means unknown and scores 1.0, matching the paper's own
convention for an unknown query rate.

**Back-pointers versus the rolling table.** The paper stores back-pointers
"alongside the DP table" while writing `dp <- new_dp`, which discards the
predecessor layer reconstruction needs. Every step table is retained instead.

### Not implemented

**Thermodynamic scoring.** The paper fixes `w_g = w_u = 0`. `deltaG` and
`deltaGerr` are parsed and stored so the scorer can be enabled without a schema
change, but no term uses them yet.

**Non-integer stoichiometry.** 308 reactions carry fractional or
scientific-notation coefficients, including biomass reactions. `stoich_matrix` is
integral, so `csv-to-json.py` skips them and reports the count.

**`filter-source` is nearly single-valued.** `grep -c KEGG` returns 0 on both
CSVs: the database has no cross-reference column, and the paper's
`KEGG:...|ModelSEED:...` strings were hand-authored for its ten-record sample.
Species filtering is the mechanism that actually reduces the database.

## Testing

```bash
cmake --build /workspaces/chemical-compiler/builds/chemcomp/build --target check-standalone
```

| Test | Covers |
|---|---|
| `db-loader-equivalence.mlir` | Both loaders produce identical databases |
| `db-slices.mlir` | Slice selection and the three filters |
| `db-errors.mlir` | Missing, malformed and inconsistent databases |
| `preprocessing-oscillator.mlir` | Paper section 5.4; the treewidth correction |
| `preprocessing-order.mlir` | Min-fill beats input order on a chain |
| `preprocessing-candidates.mlir` | Signature lookup modes; the truncation remark |
| `dp-oscillator.mlir` | Reproduces the paper's reported assignment |
| `dp-optimality.mlir` | The DP beats greedy where lookahead matters |
| `dp-injectivity.mlir` | Both injectivity modes on one input |
| `dp-determinism.mlir` | Tie-breaking and repeatability |
| `dp-no-match.mlir` | Partial coverage warning and empty-string padding |
| `scale-truncation.mlir` | 5,000-motif database with the cap engaged |
| `end-to-end-json.mlir` | Matched names reach the emulator JSON |
