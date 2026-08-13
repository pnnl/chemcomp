//===- Scoring.cpp - Match scoring ------------------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// LocalScore and its two terms, from papers/cf26-crn2db section 3.4.
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/CRNToDB/Scoring.h"

#include <cmath>

using namespace chemcomp;

double chemcomp::stoichSimilarity(const QueryCRN &query, AbstractReaction r,
                                  const Motif &motif,
                                  llvm::ArrayRef<SpeciesId> binding,
                                  llvm::ArrayRef<AbstractSpecies> support) {
  // Cosine similarity between the query row and the motif's stoichiometry
  // pulled back through the binding: S_s = (N_Q[i] . (nu_j o phi)) /
  // (||N_Q[i]|| ||nu_j o phi||).
  //
  // Both vectors are indexed over the reaction's support, so only bound
  // species contribute. A match satisfying the paper's stoichiometric
  // compatibility rule exactly scores 1.0, since the two vectors coincide.
  double dot = 0.0;
  double normQuery = 0.0;
  double normMotif = 0.0;

  for (auto [i, s] : llvm::enumerate(support)) {
    double q = query.coeff(r, s);
    double m = motif.stoichFor(binding[i]);
    dot += q * m;
    normQuery += q * q;
    normMotif += m * m;
  }

  if (normQuery <= 0.0 || normMotif <= 0.0)
    return 0.0;
  return dot / (std::sqrt(normQuery) * std::sqrt(normMotif));
}

double chemcomp::rateSimilarity(const QueryCRN &query, AbstractReaction r,
                                const Motif &motif, double sigma) {
  // S_r = exp(-(log kappa_j - log k_i)^2 / sigma^2), and 1.0 when either rate
  // is unknown -- the paper's convention for an unknown query rate, extended
  // to the motif side because the schema has no scalar rate field and a
  // missing rate_range must not be read as a disagreement.
  double queryRate = query.rateFor(r);
  double motifRate = motif.getNominalRate();

  if (std::isnan(queryRate) || std::isnan(motifRate))
    return 1.0;
  if (queryRate <= 0.0 || motifRate <= 0.0)
    return 1.0;
  if (sigma <= 0.0)
    return 1.0;

  const double delta = std::log(motifRate) - std::log(queryRate);
  return std::exp(-(delta * delta) / (sigma * sigma));
}

double chemcomp::localScore(const QueryCRN &query, AbstractReaction r,
                            const Motif &motif,
                            llvm::ArrayRef<SpeciesId> binding,
                            llvm::ArrayRef<AbstractSpecies> support,
                            const ScoringWeights &weights) {
  // LocalScore = w_c + w_s * S_s + w_r * S_r.
  //
  // w_c is a flat per-match reward, which is what makes coverage the primary
  // objective: matching one more reaction beats improving the quality of the
  // matches already made, as long as w_c dominates the other terms.
  return weights.coverage +
         weights.stoich * stoichSimilarity(query, r, motif, binding, support) +
         weights.rate * rateSimilarity(query, r, motif, weights.sigma);
}
