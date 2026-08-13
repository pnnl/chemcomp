//===- CRNMatcher.h - Abstract CRN matching ---------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// The query CRN, the preprocessing that precedes the dynamic program, and the
// match result. See papers/cf26-crn2db sections 3 and 4.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_CRNTODB_CRNMATCHER_H
#define CHEMCOMP_CONVERSION_CRNTODB_CRNMATCHER_H

#include "ChemComp/Conversion/CRNToDB/MotifDatabase.h"

#include "llvm/ADT/SmallBitVector.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/raw_ostream.h"

namespace mlir::cc {

/// Index of an abstract species (a column of stoich_matrix).
using AbstractSpecies = unsigned;
/// Index of an abstract reaction (a row of stoich_matrix).
using AbstractReaction = unsigned;

//===----------------------------------------------------------------------===//
// QueryCRN
//===----------------------------------------------------------------------===//

/// The abstract CRN to be matched: the paper's Q = (S_Q, R_Q, N_Q, k_Q).
///
/// `rows` is the stoichiometric matrix in ChemComp's convention -- reactions
/// are rows, species are columns -- read directly from the `stoich_matrix`
/// attribute of an `acr.crn_generic` operation.
struct QueryCRN {
  llvm::SmallVector<llvm::SmallVector<int32_t, 8>, 8> rows;
  unsigned numReactions = 0;
  unsigned numSpecies = 0;

  /// Optional per-reaction rate constants. Empty, or NaN at an index, means
  /// the rate is unknown and the kinetic score contributes 1.0 there.
  llvm::SmallVector<double, 8> rates;

  int32_t coeff(AbstractReaction r, AbstractSpecies s) const {
    return rows[r][s];
  }

  /// supp(r): the species with a non-zero coefficient in reaction r.
  llvm::SmallVector<AbstractSpecies, 8> support(AbstractReaction r) const;

  /// The reaction's non-zero coefficients, sorted -- its stoichiometric
  /// signature in the paper's sense.
  llvm::SmallVector<int32_t, 8> signatureCoeffs(AbstractReaction r) const;

  double rateFor(AbstractReaction r) const {
    if (r >= rates.size())
      return std::numeric_limits<double>::quiet_NaN();
    return rates[r];
  }
};

//===----------------------------------------------------------------------===//
// Preprocessing
//===----------------------------------------------------------------------===//

/// Everything computed before the dynamic program runs: the elimination
/// ordering, the per-step frontiers, and the pruned candidate domains.
struct Preprocessing {
  /// Reaction processing order, the paper's pi.
  llvm::SmallVector<AbstractReaction, 8> order;

  /// frontier[t] holds the species live *after* processing order[0..t-1], so
  /// frontier has order.size() + 1 entries and frontier[0] and the last entry
  /// are both empty. Each is sorted ascending.
  llvm::SmallVector<llvm::SmallVector<AbstractSpecies, 8>, 8> frontier;

  /// Candidate concrete species per abstract species after AC-3, sorted
  /// ascending. An empty domain means the species cannot be bound at all.
  llvm::SmallVector<llvm::SmallVector<SpeciesId, 8>, 8> domains;

  /// Candidate motifs per reaction, in ascending motif-id order.
  llvm::SmallVector<llvm::SmallVector<MotifIdx, 8>, 8> candidates;

  /// Reactions whose candidate set was truncated by the cap, with the count
  /// before truncation. Truncation makes the result heuristic rather than
  /// exact, so the caller must say so.
  llvm::SmallVector<std::pair<AbstractReaction, unsigned>, 4> truncated;

  /// True when AC-3 hit its iteration budget and stopped early, leaving
  /// domains larger than a full run would have.
  bool arcConsistencyIncomplete = false;

  /// max |F_t|, the treewidth bound W that governs the DP's state space.
  unsigned maxFrontierWidth() const;

  void print(llvm::raw_ostream &os, const QueryCRN &query,
             const MotifDatabase &db) const;
};

/// Options controlling preprocessing.
struct PreprocessOptions {
  /// Cap on candidate motifs per reaction; 0 means unlimited.
  unsigned maxCandidates = 256;

  /// Require exact signature equality instead of the coverage-correct superset
  /// lookup. Paper-literal, and much cheaper, but misses motifs carrying
  /// species beyond the query's support.
  bool strictSignature = false;

  /// Upper bound on AC-3 revision steps before giving up and keeping the
  /// domains as they stand.
  unsigned arcConsistencyBudget = 100000;
};

/// Build the species-sharing graph, derive the elimination ordering and
/// frontiers, collect candidate motifs, and prune domains with AC-3.
Preprocessing preprocess(const QueryCRN &query, const MotifDatabase &db,
                         const PreprocessOptions &options);

//===----------------------------------------------------------------------===//
// Ordering and frontier, exposed for testing
//===----------------------------------------------------------------------===//

/// Min-fill elimination order over abstract species.
///
/// Species are the variables of the interaction graph: two are adjacent when
/// some reaction's support contains both. At each step the species whose
/// elimination adds the fewest fill edges is removed, ties broken by degree
/// and then by index so the result is deterministic.
llvm::SmallVector<AbstractSpecies, 8>
computeMinFillOrder(const QueryCRN &query);

/// Reaction order derived from a species elimination order.
///
/// The paper computes min-fill over species but its dynamic program consumes
/// an ordering over reactions, without saying how one becomes the other. A
/// reaction is only fully determined once every species it touches has been
/// introduced, and that last introduction is what bounds the frontier, so
/// reactions are keyed by the latest elimination position in their support.
/// Ties break on the sorted position vector, then the reaction index.
llvm::SmallVector<AbstractReaction, 8>
deriveReactionOrder(const QueryCRN &query,
                    llvm::ArrayRef<AbstractSpecies> speciesOrder);

/// Frontiers for a reaction order, as described for Preprocessing::frontier.
llvm::SmallVector<llvm::SmallVector<AbstractSpecies, 8>, 8>
computeFrontiers(const QueryCRN &query,
                 llvm::ArrayRef<AbstractReaction> order);

} // namespace mlir::cc

#endif // CHEMCOMP_CONVERSION_CRNTODB_CRNMATCHER_H
