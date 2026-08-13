//===- FrontierDP.h - Frontier dynamic program ------------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_CRNTODB_FRONTIERDP_H
#define CHEMCOMP_CONVERSION_CRNTODB_FRONTIERDP_H

#include "ChemComp/Conversion/CRNToDB/CRNMatcher.h"
#include "ChemComp/Conversion/CRNToDB/MotifDatabase.h"
#include "ChemComp/Conversion/CRNToDB/Scoring.h"

#include <optional>

namespace chemcomp {

/// The outcome of matching one CRN.
struct MatchResult {
  /// Concrete species bound to each abstract species; nullopt where unbound.
  llvm::SmallVector<std::optional<SpeciesId>, 8> assignment;

  /// Motif matched to each reaction; nullopt where unmatched.
  llvm::SmallVector<std::optional<MotifIdx>, 8> matches;

  /// Total LocalScore over matched reactions.
  double score = 0.0;

  /// Reactions left unmatched, ascending.
  llvm::SmallVector<AbstractReaction, 4> unmatched;

  /// True when beam search discarded states, so the result may be suboptimal.
  bool beamTruncated = false;

  /// Concrete species bound to more than one abstract species.
  ///
  /// Only ever non-empty in the default frontier-local mode, where a binding
  /// is forgotten once its species leaves the frontier and can be reused. See
  /// DPOptions::globalInjectivity.
  llvm::SmallVector<SpeciesId, 4> injectivityViolations;

  unsigned numMatched() const {
    unsigned n = 0;
    for (const auto &m : matches)
      if (m)
        ++n;
    return n;
  }
};

/// Options for the dynamic program.
struct DPOptions {
  ScoringWeights weights;

  /// Beam width; 0 runs the exact dynamic program. A positive width keeps only
  /// that many states per step, ordered by the tie-breaking rule.
  unsigned beamWidth = 0;

  /// Track concrete species retired from the frontier so the final assignment
  /// is globally injective.
  ///
  /// The paper requires a globally injective phi, but its DP state cannot
  /// express that: trim() drops a binding once its species leaves the
  /// frontier, so a retired chemical can be silently reused. Enabling this
  /// carries the used set in the state key, which is correct but makes the
  /// state space exponential in the species count rather than the treewidth,
  /// destroying the complexity argument. Off by default, matching the paper.
  bool globalInjectivity = false;

  /// Cap on DP states retained per step before beam truncation applies, as a
  /// guard against pathological inputs. 0 means unlimited.
  unsigned maxStatesPerStep = 200000;
};

/// Run the frontier dynamic program of papers/cf26-crn2db Algorithm 1.
MatchResult runFrontierDP(const QueryCRN &query, const MotifDatabase &db,
                          const Preprocessing &prep, const DPOptions &options);

} // namespace chemcomp

#endif // CHEMCOMP_CONVERSION_CRNTODB_FRONTIERDP_H
