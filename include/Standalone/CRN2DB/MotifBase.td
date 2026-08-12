//===- MotifBase.td - Reaction motif database schema -------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Schema for the compiled-in reaction motif fixture consumed by the CRN2DB
// matching pass (papers/cf26-crn2db). `chemcomp-tblgen -gen-motif-db` turns
// records defined against these classes into a static C++ table.
//
// This describes the *fixture* database only -- a small curated set used by
// lit tests, so they run without any external data file. The full ~43k-motif
// ModelSEED database is loaded from JSON at runtime instead; compiling it in
// would cost minutes of build time and make the build depend on the untracked
// MotifsDB CSVs. Both paths produce the same chemcomp::MotifDatabase.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_MOTIF_BASE
#define CHEMCOMP_MOTIF_BASE

//===----------------------------------------------------------------------===//
// Motif
//===----------------------------------------------------------------------===//

// One concrete reaction from the database.
//
// `Species` and `Stoich` are positionally parallel: Stoich[i] is the net
// stoichiometric coefficient of Species[i]. Negative entries are reactants,
// positive are products, matching ChemComp's convention. Zero entries are
// dropped by the generator.
//
// TableGen has no set or map type, so the pair of parallel lists is how the
// paper's `D_j` (species set) and `nu_j : D_j -> Z` (stoichiometry map) are
// spelled here. The generated C++ materializes the real set/map semantics: a
// sorted, duplicate-free species vector queried by binary search.
//
// Real-valued fields are strings because TableGen has no float type. The
// generator parses them; the empty string means "unknown", which is a distinct
// state from any numeric value and must not collapse to 0.
class Motif<string id, list<string> sp, list<int> st> {
  // Unique identifier, e.g. "rxn00006". Sort order of the generated table.
  string MotifId = id;

  // The reaction's species, by concrete name.
  list<string> Species = sp;

  // Net stoichiometric coefficients, parallel to Species.
  list<int> Stoich = st;

  // Gibbs free energy in kJ/mol, and its uncertainty. Empty means unknown --
  // 43.1% of ModelSEED rows carry a sentinel meaning exactly that.
  //
  // Parsed and stored but unused by the matcher today: the paper fixes the
  // thermodynamic weights w_g = w_u = 0. Present so the scorer can be enabled
  // without a schema migration.
  string DeltaG = "";
  string DeltaGErr = "";

  // Plausible rate-constant bounds, as [lo, hi] decimal strings. Empty means
  // unknown. The nominal rate kappa_j is derived as the geometric mean, since
  // S_r compares rates in log space.
  list<string> RateRange = [];

  // Provenance, e.g. "ModelSEED:rxn00006" or "KEGG:R00009|ModelSEED:rxn00006".
  // Sliceable by prefix.
  string Source = "";
}

//===----------------------------------------------------------------------===//
// Extension points
//===----------------------------------------------------------------------===//

// Typed boolean properties live on subclasses rather than on Motif itself, so
// the base schema stays the paper's and new columns arrive without touching
// existing records.
//
// The generator emits every bit declared here into the flags word, so adding a
// property means adding a subclass and a line in MotifDBGen.cpp -- no change to
// MotifDatabase or the pass.
//
// SolubleMotif is the worked example establishing that pattern. Planned
// follow-ons (thermodynamic stability, further implementation-specific subsets)
// follow the same shape.
class SolubleMotif<string id, list<string> sp, list<int> st>
    : Motif<id, sp, st> {
  // True when every participant is soluble in aqueous phase, so the reaction is
  // usable by an implementation restricted to solution-phase chemistry.
  bit Soluble = 1;
}

//===----------------------------------------------------------------------===//
// Slice
//===----------------------------------------------------------------------===//

// A named subset of the database, selectable at runtime via the pass option
// `db-slice`.
//
// Slicing is not merely a test-runtime convenience. Measured over the full
// database (43,455 integer-coefficient motifs), only 913 distinct
// stoichiometric signatures exist and the single bucket (-1,-1,1,1) holds
// 11,421 motifs; 90.2% of motifs sit in buckets larger than 100. Signature
// indexing alone therefore does not bound the DP's candidate set the way the
// paper assumes, so restricting the database is what keeps matching tractable.
class Slice<string name, list<Motif> motifs> {
  string Name = name;
  list<Motif> Motifs = motifs;

  // Human-readable purpose, emitted as a comment in the generated table.
  string Summary = "";
}

#endif // CHEMCOMP_MOTIF_BASE
