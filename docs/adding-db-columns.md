# Adding a column to the motif database

This guide walks you through adding a new numeric column to the CRN2DB motif
database and using it in the matching score. It is written for someone who is
comfortable with Python and the command line but has not touched C++ or
compiler code before.

The worked example adds a **`stability`** column: a 0–100 score describing how
dynamically stable a reaction is, computed with control-theory tools (regulation,
perturbation rejection, finite-time Lyapunov exponents, ...). The same recipe
works for any other real-valued per-reaction column.

Background on the pass itself is in [crn2db.md](crn2db.md).

## 0. What you are changing

A column travels through five layers. Each layer is one section below.

```text
 stability CSV ──► utils/*.py ──► JSON file ─────────────┐ (runtime loader)
                       │                                  ▼
                       └──► MotifsFixture.td ──► chemcomp-tblgen ──► MotifsFixture.inc
                                                                           │ (compiled-in loader)
                                                                           ▼
                                              Motif struct (MotifDatabase) ──► Scoring.cpp ──► -crn-db-match
```

There are two ways a database reaches the compiler, and both must learn about
the new column, otherwise they stop being equivalent:

| Path | Used when | Data comes from |
|---|---|---|
| Compiled-in fixture | `-crn-db-match='db-slice=...'` | [MotifsFixture.td](../include/ChemComp/Conversion/CRNToDB/MotifsFixture.td), turned into C++ by `chemcomp-tblgen` |
| JSON file | `-crn-db-match='db=path.json'` | A file you generate with the scripts in [utils/](../utils/) |

Files you will touch (checklist at the end):

| Layer | Files |
|---|---|
| Data scripts | [motifsdb_parse.py](../utils/motifsdb_parse.py), [csv-to-json.py](../utils/csv-to-json.py), [csv-to-td.py](../utils/csv-to-td.py) |
| TableGen schema | [MotifBase.td](../include/ChemComp/Conversion/CRNToDB/MotifBase.td), [MotifsFixture.td](../include/ChemComp/Conversion/CRNToDB/MotifsFixture.td) |
| TableGen backend | [MotifDBGen.cpp](../tools/chemcomp-tblgen/MotifDBGen.cpp) |
| In-memory row + loaders | [MotifDatabase.h](../include/ChemComp/Conversion/CRNToDB/MotifDatabase.h), [MotifDatabase.cpp](../lib/Conversion/CRNToDB/MotifDatabase.cpp), [MotifDatabaseJSON.cpp](../lib/Conversion/CRNToDB/MotifDatabaseJSON.cpp), [MotifDatabaseFixture.cpp](../lib/Conversion/CRNToDB/MotifDatabaseFixture.cpp) |
| Scoring | [Scoring.h](../include/ChemComp/Conversion/CRNToDB/Scoring.h), [Scoring.cpp](../lib/Conversion/CRNToDB/Scoring.cpp), [Passes.td](../include/ChemComp/Conversion/Passes.td), [CRNDBMatchPass.cpp](../lib/Conversion/CRNToDB/CRNDBMatchPass.cpp) |
| Tests | [tiny-db.json](../test/Conversion/CRNToDB/Inputs/tiny-db.json), [db-loader-equivalence.mlir](../test/Conversion/CRNToDB/db-loader-equivalence.mlir), [end-to-end-json.mlir](../test/Conversion/CRNToDB/end-to-end-json.mlir) |

Every command below is run from the `chemcomp/` directory (the repository root
of the public project).

## 1. Set up your machine (macOS)

1. Install [Pixi](https://pixi.sh):

   ```sh
   curl -fsSL https://pixi.sh/install.sh | bash
   ```

2. Pixi currently lists only `linux-64` under `platforms` in
   [pixi.toml](../pixi.toml). Add your Mac's platform (`osx-arm64` for Apple
   silicon, `osx-64` for Intel):

   ```toml
   platforms = ["linux-64", "osx-arm64"]
   ```

   > **Untested.** macOS is not yet part of the project's CI. Two things may
   > need adjusting: the conda-forge LLVM/MLIR 18.1.2 package must exist for your
   > platform, and the `link-llvm-utils` task in `pixi.toml` links `FileCheck`
   > and friends from `$CONDA_PREFIX/libexec/llvm`, which may differ on macOS.
   > If either fails, ask for help before changing anything else.

3. Install and build:

   ```sh
   pixi install          # re-solves and updates pixi.lock
   pixi run build        # configures build/pixi-<env> and builds chemcomp-opt
   pixi run check        # builds everything and runs the whole test-suite
   ```

**Expected result:** `pixi run check` ends with all tests passing. Do not
continue until this baseline is green; otherwise you cannot tell your own
breakage from pre-existing breakage.

> `pixi run build` only builds `chemcomp-opt`. The tests also need
> `chemcomp-translate` and `chemcomp-tblgen`, which `pixi run check` builds.

## 2. Get your data

The raw reaction database is a pair of CSV files, `reactions_reduced.csv`
(`motif_id,equation,deltag,deltagerr`) and `compounds_reduced.csv`, in the
`MotifsDB` submodule at the repository root. The stability extension is a
separate CSV keyed by the same `motif_id`:

```text
MotifsDB/stability_extension/canonical_dynamic_reliability_10motifs_equal_weights_v0_1.csv
```

It has one row per motif and, among others, the columns `s_reg`, `s_pert`,
`s_bif` (often `N/A`), `s_FTLE`, `s_num`, and `equal_weight_scalar_score`
(0–100, the equal-weighted mean of the available components). We use
`equal_weight_scalar_score` as the `stability` column. Only 10 motifs
(`rxn00001`–`rxn00011`) have a score; all others must remain *unknown*, which the
database represents as a missing value (not `0`).

**Producing your own scores.** The reservoir emulator in
[chemcomp-reservoir-emulator](../tools/chemcomp-reservoir-emulator/) can simulate a
network such as
[oscillatory.json](../tools/chemcomp-reservoir-emulator/network_configs/oscillatory.json)
(`pixi run emulator --network network_configs/oscillatory.json ...`), and its
`chemical_reservoir/metrics.py` and `chemical_reservoir/sensitivity.py` compute
error-growth and sensitivity numbers. It does **not** compute Jacobians,
eigenvalues or Lyapunov exponents, so a control-theory score needs its own
script. Whatever you write only has to emit a CSV with a `motif_id` column and one
numeric column.

## 3. Step 1 — Python: carry the column through the converters

The scripts in [utils/](../utils/) convert the CSV database into the two formats
the compiler reads. They share one parser,
[motifsdb_parse.py](../utils/motifsdb_parse.py).

**3a. Parser.** Add an optional field to the `Motif` dataclass and a merge helper:

```python
@dataclass
class Motif:
    motif_id: str
    species: List[str]
    stoich: List[int]
    delta_g: Optional[float]
    delta_g_err: Optional[float]
    source: str
    stability: Optional[float] = None   # new


def load_stability(path: str, column: str = "equal_weight_scalar_score") -> dict:
    """motif_id -> score; 'N/A' and empty cells are omitted (= unknown)."""
    import csv
    scores = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            try:
                scores[row["motif_id"]] = float(row[column])
            except ValueError:
                pass
    return scores
```

**3b. JSON export.** In [csv-to-json.py](../utils/csv-to-json.py) add a
`--stability-csv` argument, and write the key only when the value is known (next
to where `deltaG` is written):

```python
ap.add_argument("--stability-csv", help="CSV with motif_id and a score column")
...
if m.stability is not None:
    record["stability"] = m.stability
```

(After `parse_database`, set `m.stability = scores.get(m.motif_id)` for each
motif.)

**3c. TableGen export.** In `emit_motif` of [csv-to-td.py](../utils/csv-to-td.py):

```python
if motif.stability is not None:
    print(f'  let Stability = "{motif.stability}";', file=out)
```

**Run it:**

```sh
./utils/csv-to-json.py \
    --reactions MotifsDB/reactions_reduced.csv \
    --compounds MotifsDB/compounds_reduced.csv \
    --stability-csv MotifsDB/stability_extension/canonical_dynamic_reliability_10motifs_equal_weights_v0_1.csv \
    -o motifs.json
grep -m1 '"stability"' motifs.json

./utils/csv-to-td.py \
    --reactions MotifsDB/reactions_reduced.csv \
    --compounds MotifsDB/compounds_reduced.csv \
    --stability-csv MotifsDB/stability_extension/canonical_dynamic_reliability_10motifs_equal_weights_v0_1.csv \
    --id rxn00001 --id rxn00006 --slice-name StabilitySlice
```

**Expected result:** `rxn00001` has `"stability": 88.654242`; the `.td` output
contains `let Stability = "88.654242";`. Motifs without a score have no
`stability` key at all. (Pass `--stability-csv` to `csv-to-td.py` too, mirroring
the JSON script.)

## 4. Step 2 — TableGen: declare the column

TableGen is LLVM's data-description language. The database schema is a TableGen
`class`; each reaction is a `def` of that class.

In [MotifBase.td](../include/ChemComp/Conversion/CRNToDB/MotifBase.td), add a
field to `class Motif` next to `DeltaG`:

```tablegen
  // Dynamic-stability score in [0, 100]. Empty means unknown.
  string Stability = "";
```

It is a `string` because TableGen has no floating point type; the backend parses
it. The empty string means unknown.

In [MotifsFixture.td](../include/ChemComp/Conversion/CRNToDB/MotifsFixture.td),
set it on the motifs that have data:

```tablegen
def Rxn00001 : Motif<"rxn00001",
    ["H2O", "PPi", "Phosphate", "H+"], [-1, -1, 2, 1]> {
  let DeltaG = "-3.46";
  let DeltaGErr = "0.05";
  let RateRange = ["0.01", "1.0"];
  let Stability = "88.654242";
  let Source = "KEGG:R00004|ModelSEED:rxn00001";
}
```

The default (`""`) keeps every other record valid unchanged. Do **not** add a new
`Slice` yet: [db-errors.mlir](../test/Conversion/CRNToDB/db-errors.mlir) hard-codes
the list of slice names and would fail.

## 5. Step 3 — TableGen backend: generate the C++ row

[MotifDBGen.cpp](../tools/chemcomp-tblgen/MotifDBGen.cpp) reads the records and
prints a C++ table, `MotifsFixture.inc`. Two edits, and **their order must
match**, because rows are positional initialisers.

1. In the emitted `struct FixtureMotif`, add the field just before `flags`:

   ```cpp
   os << "  const char *source;\n";
   os << "  double stability;\n";      // new
   os << "  unsigned flags;\n";
   ```

2. In `emitMotif`, emit the value after the source and before the flags:

   ```cpp
   emitStringLiteral(os, rec->getValueAsString("Source"));
   os << ", ";
   emitOptionalDouble(os, rec->getValueAsString("Stability"), rec, "Stability");
   // ... existing code that prints the flags follows
   ```

   `emitOptionalDouble` already prints `kUnknown` (NaN) for `""` and reports a
   readable error for text that is not a number.

**Run it** (the first command regenerates the `.inc` as part of the build):

```sh
pixi run build
find build -name MotifsFixture.inc | head -1      # path of the generated file
grep -n "88.654242" "$(find build -name MotifsFixture.inc | head -1)"
```

**Expected result:** the build fails with a compile error in
`MotifDatabaseFixture.cpp` only if you reorder the struct; otherwise it succeeds
and the grep shows `rxn00001`'s row. Then update the tablegen unit tests: their
`CHECK` lines contain each row's trailing fields.

```sh
pixi run -e compiler lit -v build/pixi-compiler/test/chemcomp-tblgen
```

Fix the expectations in
[motif-db-gen.td](../test/chemcomp-tblgen/motif-db-gen.td) and
[fixture-db.td](../test/chemcomp-tblgen/fixture-db.td) by inserting the new
value (`nan` for unknown, the number otherwise) in the position you chose.

> The `lit` invocation above depends on your environment name; if the path
> doesn't exist, run `ls build` and use the `build/pixi-*` directory you see.

## 6. Step 4 — C++: store the value on a motif

This is the first real C++ edit but it is mechanical: add a field, thread it
through two functions.

**6a. The row.** In
[MotifDatabase.h](../include/ChemComp/Conversion/CRNToDB/MotifDatabase.h), next to
`deltaGErr`:

```cpp
/// Dynamic-stability score in [0, 100]; NaN when unknown.
double stability = std::numeric_limits<double>::quiet_NaN();
```

and add a parameter at the end of `addMotif`:

```cpp
bool addMotif(llvm::StringRef id, llvm::ArrayRef<llvm::StringRef> speciesNames,
              llvm::ArrayRef<int32_t> stoich, double deltaG, double deltaGErr,
              llvm::ArrayRef<double> rateRange, llvm::StringRef source,
              double stability = std::numeric_limits<double>::quiet_NaN());
```

The default keeps other callers compiling. Mirror the parameter in
[MotifDatabase.cpp](../lib/Conversion/CRNToDB/MotifDatabase.cpp) (`addMotif`,
around line 108: store it in the new `Motif`).

**6b. The dump.** In `MotifDatabase::dump` of the same file, append the value at
the **end** of the motif line, so the existing `CHECK-NEXT` lines in other tests
(which match a prefix of the line) keep passing:

```cpp
if (std::isnan(m.stability))
  os << " stability=unknown";
else
  os << llvm::format(" stability=%g", m.stability);
```

(Place it after the existing `source=` print; read the end of the loop body to
find the exact spot.)

**6c. The two loaders.**

- JSON, [MotifDatabaseJSON.cpp](../lib/Conversion/CRNToDB/MotifDatabaseJSON.cpp):

  ```cpp
  double stability;
  if (llvm::Error e = readOptionalDouble(*obj, "stability", index, stability))
    return e;
  ...
  db.addMotif(*id, speciesNames, stoich, deltaG, deltaGErr, rateRange, source,
              stability)
  ```

  `readOptionalDouble` already treats a missing key and `null` as unknown, and
  rejects a non-numeric value with an error that names the motif.

- Fixture, [MotifDatabaseFixture.cpp](../lib/Conversion/CRNToDB/MotifDatabaseFixture.cpp):
  pass `row.stability` as the last argument of `db.addMotif(...)`.

**Run it:**

```sh
pixi run build
grep -rn "addMotif(" lib include tools test     # every caller should be updated or use the default
```

## 7. Step 5 — Scoring: use the column

The score of one match is computed in
[Scoring.cpp](../lib/Conversion/CRNToDB/Scoring.cpp):
`LocalScore = w_c + w_s·S_s + w_r·S_r` (coverage, stoichiometric similarity,
kinetic proximity). You add a fourth term `w_t·S_t`, where `S_t` is the motif's
stability rescaled to `[0, 1]`.

**7a. Weight.** In [Scoring.h](../include/ChemComp/Conversion/CRNToDB/Scoring.h):

```cpp
struct ScoringWeights {
  double coverage = 2.0; // w_c
  double stoich = 1.0;   // w_s
  double rate = 0.5;     // w_r
  double sigma = 1.0;    // width of the log-space Gaussian in S_r
  double stability = 0.0; // w_t; 0 keeps every existing score unchanged
};

/// S_t: stability / 100, and 1.0 when the motif's stability is unknown.
double stabilitySimilarity(const Motif &motif);
```

Keep the default weight at `0.0`. This is important: existing tests then give
bit-identical results, and your change only has an effect when asked for.

**7b. Term.** In [Scoring.cpp](../lib/Conversion/CRNToDB/Scoring.cpp):

```cpp
double mlir::cc::stabilitySimilarity(const Motif &motif) {
  // Unknown scores are neutral (1.0), like an unknown rate in rateSimilarity:
  // a missing measurement must not be read as instability.
  if (std::isnan(motif.stability))
    return 1.0;
  return std::min(std::max(motif.stability / 100.0, 0.0), 1.0);
}
```

(add `#include <algorithm>`), and extend `localScore`:

```cpp
return weights.coverage +
       weights.stoich * stoichSimilarity(query, r, motif, binding, support) +
       weights.rate * rateSimilarity(query, r, motif, weights.sigma) +
       weights.stability * stabilitySimilarity(motif);
```

> **Design choice.** Treating unknown as `1.0` means motifs *without* a score
> are never penalised. If you would rather prefer measured-stable motifs over
> unmeasured ones, return a lower neutral value (e.g. `0.5`) instead. Pick one
> and document it in the comment.

**7c. Command-line option.** In [Passes.td](../include/ChemComp/Conversion/Passes.td),
beside `wr`:

```tablegen
Option<"weightStability", "wt", "double", /*default=*/"0.0",
       "Weight of the dynamic-stability term">,
```

and in [CRNDBMatchPass.cpp](../lib/Conversion/CRNToDB/CRNDBMatchPass.cpp), beside
the other weights (around line 153):

```cpp
dpOptions.weights.stability = weightStability;
```

The scorer is called from
[FrontierDP.cpp](../lib/Conversion/CRNToDB/FrontierDP.cpp) via `localScore`, which
already receives the whole `ScoringWeights`, so nothing else changes.

**Run it:**

```sh
pixi run build
# default behaviour: unchanged
pixi run chemcomp-opt -crn-db-match='db-slice=OscillatorSlice' \
    test/Conversion/CRNToDB/end-to-end-json.mlir
# stability enabled: same input, different weight
pixi run chemcomp-opt -crn-db-match='db-slice=OscillatorSlice wt=1.0' \
    test/Conversion/CRNToDB/end-to-end-json.mlir
```

## 8. Step 6 — Test

Tests are `lit` files: an MLIR input with `RUN:` lines (commands) and `CHECK:`
lines (expected output, matched by `FileCheck`).

**8a. Loader equivalence.** Add `"stability"` to the records in
[tiny-db.json](../test/Conversion/CRNToDB/Inputs/tiny-db.json), set matching
`let Stability = ...` on the `TinySlice` records in `MotifsFixture.td`, and extend
the `CHECK-NEXT` lines of
[db-loader-equivalence.mlir](../test/Conversion/CRNToDB/db-loader-equivalence.mlir):

```mlir
// CHECK-NEXT: motif rxn00001 ... source=KEGG:R00004|ModelSEED:rxn00001 stability=88.6542
```

Because the same `CHECK` is run against the fixture *and* the JSON file, this one
test proves both loaders agree.

```sh
pixi run chemcomp-opt \
    -crn-db-match='db=test/Conversion/CRNToDB/Inputs/tiny-db.json dump-db=1' \
    test/Conversion/CRNToDB/db-loader-equivalence.mlir
```

**8b. A scoring test.** Copy the shape of
[end-to-end-json.mlir](../test/Conversion/CRNToDB/end-to-end-json.mlir): one
input, two `RUN:` lines, one with `wt=0` and one with a large `wt`, using a
database JSON in `Inputs/` where two candidate motifs have the same stoichiometry
but different stability. Check that the more stable motif's species names appear
in the output only when `wt` is large.

**8c. Run the tests.**

```sh
# one test (path depends on the environment name under build/)
pixi run -e compiler lit -v build/pixi-compiler/test/Conversion/CRNToDB/db-loader-equivalence.mlir
# the whole CRNToDB directory
pixi run -e compiler lit -v build/pixi-compiler/test/Conversion/CRNToDB
# everything
pixi run check
```

**Expected result:** all pass. `lit` prints `Passed: N`; any `FAIL` shows the
actual output next to the `CHECK` that did not match.

**8d. Lint and format.**

```sh
pixi run fmt-py && pixi run lint-py        # Python
pixi run fmt-cpp && pixi run lint-cpp      # C++
```

## 9. Troubleshooting

| Symptom | Likely cause |
|---|---|
| Compile error in `MotifsFixture.inc` about too many/few initialisers | Struct field order in `MotifDBGen.cpp` does not match the order `emitMotif` prints values |
| `chemcomp-tblgen: ... Stability ...` error | A `let Stability = "..."` is not a plain decimal number |
| `db-errors.mlir` fails with a different slice list | A new `Slice` was added; update the list on that test's `CHECK` line |
| `deltaG` and others in old tests now fail | The dump line was changed in the middle; append new fields at the end |
| JSON loads but `stability=unknown` everywhere | Key spelled differently, or the value is `1e7` or larger (the loader treats ≥ 1e7 as the "unknown" sentinel) |
| Scores unchanged with `wt=1` | `weightStability` not copied in `CRNDBMatchPass.cpp`, or the candidate motifs have no stability (neutral term) |
| `lit` cannot find the test | Wrong `build/pixi-*` directory; run `ls build` |

## 10. Checklist

- [ ] `pixi run check` green before you start
- [ ] [motifsdb_parse.py](../utils/motifsdb_parse.py), [csv-to-json.py](../utils/csv-to-json.py), [csv-to-td.py](../utils/csv-to-td.py)
- [ ] [MotifBase.td](../include/ChemComp/Conversion/CRNToDB/MotifBase.td), [MotifsFixture.td](../include/ChemComp/Conversion/CRNToDB/MotifsFixture.td)
- [ ] [MotifDBGen.cpp](../tools/chemcomp-tblgen/MotifDBGen.cpp) (struct **and** emit, same order)
- [ ] [MotifDatabase.h](../include/ChemComp/Conversion/CRNToDB/MotifDatabase.h) / [.cpp](../lib/Conversion/CRNToDB/MotifDatabase.cpp) (field, `addMotif`, `dump` at end of line)
- [ ] [MotifDatabaseJSON.cpp](../lib/Conversion/CRNToDB/MotifDatabaseJSON.cpp), [MotifDatabaseFixture.cpp](../lib/Conversion/CRNToDB/MotifDatabaseFixture.cpp)
- [ ] [Scoring.h](../include/ChemComp/Conversion/CRNToDB/Scoring.h), [Scoring.cpp](../lib/Conversion/CRNToDB/Scoring.cpp), default weight `0.0`
- [ ] [Passes.td](../include/ChemComp/Conversion/Passes.td), [CRNDBMatchPass.cpp](../lib/Conversion/CRNToDB/CRNDBMatchPass.cpp)
- [ ] Tests updated and added; `pixi run check`, `lint-py`, `lint-cpp` green
