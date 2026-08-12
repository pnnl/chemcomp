//===- CRN2DBPasses.td - CRN2DB passes ---------------------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CRN2DB_PASSES
#define CHEMCOMP_CRN2DB_PASSES

include "mlir/Pass/PassBase.td"

def CRNDBMatch : Pass<"crn-db-match", "::mlir::ModuleOp"> {
  let summary = "Match abstract CRNs against a concrete reaction database";
  let description = [{
    Binds the anonymous species of an `acr.crn_generic` operation to concrete
    chemicals, and each of its reactions to a database reaction motif with
    compatible stoichiometry.

    The pass reads the `stoich_matrix` attribute, matches it against a motif
    database, and annotates the operation in place with two discardable
    attributes:

      * `species_names` -- concrete name per abstract species column, and
      * `motif_ids` -- matched motif id per reaction row.

    Unmatched entries are the empty string. When coverage is partial the pass
    emits a warning naming the unmatched reactions and proceeds with the best
    partial matching, so downstream passes can still consume the result.

    The database comes either from the compiled-in TableGen fixture (the
    default, used by tests so they need no external file) or from a JSON file
    given by `db`.

    Example:

    ```mlir
    // RUN: standalone-opt -crn-db-match='db-slice=OscillatorSlice'
    %out = acr.crn_generic {stoich_matrix = dense<...> : tensor<4x6xi32>}
        (%species, %rates) : (...) -> ...
    ```
  }];

  let options = [
    Option<"dbPath", "db", "std::string", /*default=*/"\"\"",
           "Path to a JSON motif database; empty selects the compiled-in "
           "fixture">,
    Option<"dbSlice", "db-slice", "std::string", /*default=*/"\"\"",
           "Named slice of the compiled-in fixture; empty selects every "
           "fixture motif">,
    ListOption<"filterSpecies", "filter-species", "std::string",
           "Keep only motifs containing at least one of these species">,
    Option<"filterSource", "filter-source", "std::string", /*default=*/"\"\"",
           "Keep only motifs whose source starts with this prefix">,
    Option<"maxCandidates", "max-candidates", "unsigned", /*default=*/"256",
           "Cap on candidate motifs considered per reaction (0 = unlimited). "
           "Truncation keeps the lowest motif ids and emits a remark, because "
           "it makes the result heuristic rather than exact">,
    Option<"maxPerSignature", "max-per-signature", "unsigned",
           /*default=*/"0",
           "Cap on motifs loaded per stoichiometric signature (0 = unlimited). "
           "Unlike 'max-candidates' this shrinks the database itself, so the "
           "discarded motifs are never considered and no remark is emitted">,
    Option<"weightCoverage", "wc", "double", /*default=*/"2.0",
           "Weight of the per-match coverage reward">,
    Option<"weightStoich", "ws", "double", /*default=*/"1.0",
           "Weight of the stoichiometric similarity term">,
    Option<"weightRate", "wr", "double", /*default=*/"0.5",
           "Weight of the kinetic proximity term">,
    Option<"rateSigma", "sigma", "double", /*default=*/"1.0",
           "Width of the log-space Gaussian used by the kinetic term">,
    Option<"beamWidth", "beam-width", "unsigned", /*default=*/"0",
           "Beam width; 0 runs the exact dynamic program">,
    Option<"strictSignature", "strict-signature", "bool", /*default=*/"false",
           "Require exact signature equality instead of the coverage-correct "
           "superset lookup">,
    Option<"globalInjectivity", "global-injectivity", "bool",
           /*default=*/"false",
           "Track retired species bindings so the assignment is globally "
           "injective, at the cost of a state space exponential in the species "
           "count rather than the treewidth">,
    Option<"dumpDatabase", "dump-db", "bool", /*default=*/"false",
           "Print the loaded database and stop; for debugging and loader "
           "tests">,
    Option<"dumpPreprocessing", "dump-preprocessing", "bool",
           /*default=*/"false",
           "Print the elimination ordering, frontiers, candidate sets and "
           "pruned domains instead of matching">,
  ];
}

#endif // CHEMCOMP_CRN2DB_PASSES
