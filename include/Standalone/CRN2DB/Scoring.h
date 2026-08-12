//===- Scoring.h - Match scoring --------------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef STANDALONE_CRN2DB_SCORING_H
#define STANDALONE_CRN2DB_SCORING_H

#include "Standalone/CRN2DB/CRNMatcher.h"
#include "Standalone/CRN2DB/MotifDatabase.h"

#include "llvm/ADT/ArrayRef.h"

namespace chemcomp {

/// Weights of the objective in papers/cf26-crn2db section 3.4.
///
/// The paper's defaults make coverage dominate: matching one more reaction is
/// worth more than improving matches already made.
///
/// The thermodynamic weights w_g and w_u are fixed at 0 there and are not
/// implemented; deltaG is parsed and stored so they can be enabled without a
/// schema change.
struct ScoringWeights {
  double coverage = 2.0; // w_c
  double stoich = 1.0;   // w_s
  double rate = 0.5;     // w_r
  double sigma = 1.0;    // width of the log-space Gaussian in S_r
};

/// S_s: cosine similarity between the query row and the motif's stoichiometry
/// pulled back through the binding. `binding[i]` is the concrete species bound
/// to `support[i]`.
double stoichSimilarity(const QueryCRN &query, AbstractReaction r,
                        const Motif &motif, llvm::ArrayRef<SpeciesId> binding,
                        llvm::ArrayRef<AbstractSpecies> support);

/// S_r: log-space Gaussian on the rate constants, 1.0 when either side's rate
/// is unknown.
double rateSimilarity(const QueryCRN &query, AbstractReaction r,
                      const Motif &motif, double sigma);

/// LocalScore = w_c + w_s * S_s + w_r * S_r.
double localScore(const QueryCRN &query, AbstractReaction r, const Motif &motif,
                  llvm::ArrayRef<SpeciesId> binding,
                  llvm::ArrayRef<AbstractSpecies> support,
                  const ScoringWeights &weights);

} // namespace chemcomp

#endif // STANDALONE_CRN2DB_SCORING_H
