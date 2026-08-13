//===- FrontierDP.cpp - Frontier dynamic program ----------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Algorithm 1 of papers/cf26-crn2db.
//
// Two implementation points the paper leaves open:
//
//  * It stores back-pointers "alongside the DP table" but writes the table as
//    a rolling `dp <- new_dp`, which discards the predecessor layer. Here every
//    step table is retained, so reconstruction is a walk backwards. Memory is
//    O(m * states), which is nothing at the sizes involved.
//
//  * It specifies no tie-breaking, yet results are compared by FileCheck. A
//    total order is imposed: score, then match count, then the motif-id
//    sequence, then the binding. Candidates are always visited in ascending
//    motif-id order and every domain is a sorted vector, so no iteration
//    depends on hash order.
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/CRNToDB/FrontierDP.h"

#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/MapVector.h"
#include "llvm/ADT/STLExtras.h"

#include <algorithm>

using namespace chemcomp;

namespace {

/// Scores within this distance compare equal, deferring to the tie-break.
constexpr double kScoreEpsilon = 1e-9;

/// One binding in a DP state: an abstract species and the chemical bound to it.
struct Binding {
  AbstractSpecies species;
  SpeciesId concrete;

  bool operator==(const Binding &o) const {
    return species == o.species && concrete == o.concrete;
  }
  bool operator<(const Binding &o) const {
    if (species != o.species)
      return species < o.species;
    return concrete < o.concrete;
  }
};

/// A DP state: the frontier assignment phi_t, plus the retired chemicals when
/// global injectivity is enforced.
struct StateKey {
  llvm::SmallVector<Binding, 4> bindings; // sorted by species
  llvm::SmallVector<SpeciesId, 4> retired; // sorted; empty unless enforcing

  bool operator==(const StateKey &o) const {
    return bindings == o.bindings && retired == o.retired;
  }
};

struct StateKeyInfo {
  static StateKey getEmptyKey() {
    StateKey k;
    k.bindings.push_back({~0u, ~0u});
    return k;
  }
  static StateKey getTombstoneKey() {
    StateKey k;
    k.bindings.push_back({~0u - 1, ~0u - 1});
    return k;
  }
  static unsigned getHashValue(const StateKey &k) {
    llvm::hash_code h = llvm::hash_value(k.bindings.size());
    for (const Binding &b : k.bindings)
      h = llvm::hash_combine(h, b.species, b.concrete);
    for (SpeciesId s : k.retired)
      h = llvm::hash_combine(h, s);
    return static_cast<unsigned>(h);
  }
  static bool isEqual(const StateKey &a, const StateKey &b) { return a == b; }
};

/// A reachable state and how it was reached.
struct StateEntry {
  double score = 0.0;
  unsigned numMatched = 0;
  /// Index into the previous step's table; ~0u at step 0.
  unsigned previous = ~0u;
  /// Motif matched at this step, or nullopt when the reaction was skipped.
  std::optional<MotifIdx> motif;
  /// Bindings introduced at this step, needed to rebuild the full assignment
  /// after trim() has dropped them from the state.
  llvm::SmallVector<Binding, 4> introduced;
  /// Motif ids along the path, for deterministic tie-breaking.
  llvm::SmallVector<llvm::StringRef, 4> motifPath;
};

using StepTable = llvm::MapVector<StateKey, StateEntry, // deterministic order
                                  llvm::DenseMap<StateKey, unsigned, StateKeyInfo>,
                                  llvm::SmallVector<std::pair<StateKey, StateEntry>, 8>>;

/// Total order making the DP deterministic where scores tie.
///
/// Prefers, in order: higher score, more matched reactions, the
/// lexicographically smaller motif-id sequence, then the smaller binding.
bool isBetter(const StateEntry &a, const StateEntry &b) {
  if (std::abs(a.score - b.score) > kScoreEpsilon)
    return a.score > b.score;
  if (a.numMatched != b.numMatched)
    return a.numMatched > b.numMatched;
  if (a.motifPath.size() != b.motifPath.size())
    return a.motifPath.size() < b.motifPath.size();
  for (auto [x, y] : llvm::zip(a.motifPath, b.motifPath))
    if (x != y)
      return x < y;
  return false;
}

void updateState(StepTable &table, const StateKey &key, StateEntry entry) {
  auto it = table.find(key);
  if (it == table.end()) {
    table.insert({key, std::move(entry)});
    return;
  }
  if (isBetter(entry, it->second))
    it->second = std::move(entry);
}

/// Project a state onto the frontier, retiring what leaves it.
StateKey trim(const StateKey &state, llvm::ArrayRef<AbstractSpecies> frontier,
              bool globalInjectivity) {
  StateKey out;
  out.retired = state.retired;

  for (const Binding &b : state.bindings) {
    if (std::binary_search(frontier.begin(), frontier.end(), b.species)) {
      out.bindings.push_back(b);
    } else if (globalInjectivity) {
      // The binding leaves the state but its chemical stays spoken for.
      out.retired.push_back(b.concrete);
    }
  }

  if (globalInjectivity) {
    llvm::sort(out.retired);
    out.retired.erase(std::unique(out.retired.begin(), out.retired.end()),
                      out.retired.end());
  }
  return out;
}

/// True when `concrete` is already used by this state.
bool isTaken(const StateKey &state, SpeciesId concrete) {
  for (const Binding &b : state.bindings)
    if (b.concrete == concrete)
      return true;
  return std::binary_search(state.retired.begin(), state.retired.end(),
                            concrete);
}

/// Enumerate injective extensions binding a reaction's free species to motif
/// species with the required coefficients.
///
/// `binding` is filled positionally against `support`. Bound species are fixed
/// and only checked; free species are enumerated over the motif's own species
/// carrying the right coefficient, which is what keeps this O(|Sigma_l|^d)
/// rather than O(|Sigma|^d).
void enumerateExtensions(
    const QueryCRN &query, AbstractReaction reaction, const Motif &motif,
    llvm::ArrayRef<AbstractSpecies> support, const StateKey &state,
    const Preprocessing &prep, size_t position,
    llvm::SmallVectorImpl<SpeciesId> &binding,
    llvm::SmallVectorImpl<Binding> &introduced,
    llvm::function_ref<void(llvm::ArrayRef<SpeciesId>,
                            llvm::ArrayRef<Binding>)>
        emit) {
  if (position == support.size()) {
    emit(binding, introduced);
    return;
  }

  const AbstractSpecies species = support[position];
  const int32_t want = query.coeff(reaction, species);

  // Already bound by an earlier reaction: the choice is fixed, so only check
  // that the motif offers it with the right coefficient.
  for (const Binding &b : state.bindings) {
    if (b.species != species)
      continue;
    if (motif.stoichFor(b.concrete) != want)
      return;
    binding.push_back(b.concrete);
    enumerateExtensions(query, reaction, motif, support, state, prep,
                        position + 1, binding, introduced, emit);
    binding.pop_back();
    return;
  }

  // Free: try each motif species with the required coefficient, in the motif's
  // sorted species order so enumeration is deterministic.
  for (auto [i, concrete] : llvm::enumerate(motif.species)) {
    if (motif.stoich[i] != want)
      continue;

    // The AC-3 domain is the pruned set of chemicals this species can take.
    const auto &domain = prep.domains[species];
    if (!std::binary_search(domain.begin(), domain.end(), concrete))
      continue;

    // Injectivity, within the state and within this extension.
    if (isTaken(state, concrete))
      continue;
    if (llvm::is_contained(binding, concrete))
      continue;

    binding.push_back(concrete);
    introduced.push_back({species, concrete});
    enumerateExtensions(query, reaction, motif, support, state, prep,
                        position + 1, binding, introduced, emit);
    introduced.pop_back();
    binding.pop_back();
  }
}

} // namespace

MatchResult chemcomp::runFrontierDP(const QueryCRN &query,
                                    const MotifDatabase &db,
                                    const Preprocessing &prep,
                                    const DPOptions &options) {
  const unsigned m = static_cast<unsigned>(prep.order.size());

  MatchResult result;
  result.assignment.assign(query.numSpecies, std::nullopt);
  result.matches.assign(query.numReactions, std::nullopt);

  // Every step table is kept so back-pointers stay valid; the paper's rolling
  // table cannot support the reconstruction it calls for.
  llvm::SmallVector<StepTable, 8> tables;
  tables.resize(m + 1);

  // dp[0][{}] = 0.
  tables[0].insert({StateKey{}, StateEntry{}});

  for (unsigned t = 0; t < m; ++t) {
    const AbstractReaction reaction = prep.order[t];
    const llvm::SmallVector<AbstractSpecies, 8> support =
        query.support(reaction);
    llvm::ArrayRef<AbstractSpecies> frontier = prep.frontier[t + 1];
    StepTable &next = tables[t + 1];

    for (auto [stateIndex, item] : llvm::enumerate(tables[t])) {
      const StateKey &state = item.first;
      const StateEntry &entry = item.second;

      // Option A: skip the reaction. Adds no score and no penalty.
      {
        StateEntry skipped;
        skipped.score = entry.score;
        skipped.numMatched = entry.numMatched;
        skipped.previous = static_cast<unsigned>(stateIndex);
        skipped.motif = std::nullopt;
        skipped.motifPath = entry.motifPath;
        updateState(next, trim(state, frontier, options.globalInjectivity),
                    std::move(skipped));
      }

      // Option B: match the reaction to a candidate motif.
      for (MotifIdx candidate : prep.candidates[reaction]) {
        const Motif &motif = db.getMotif(candidate);

        llvm::SmallVector<SpeciesId, 8> binding;
        llvm::SmallVector<Binding, 4> introduced;
        enumerateExtensions(
            query, reaction, motif, support, state, prep, 0, binding, introduced,
            [&](llvm::ArrayRef<SpeciesId> bound,
                llvm::ArrayRef<Binding> added) {
              StateKey extended = state;
              extended.bindings.append(added.begin(), added.end());
              llvm::sort(extended.bindings);

              StateEntry matched;
              matched.score =
                  entry.score + localScore(query, reaction, motif, bound,
                                           support, options.weights);
              matched.numMatched = entry.numMatched + 1;
              matched.previous = static_cast<unsigned>(stateIndex);
              matched.motif = candidate;
              matched.introduced.assign(added.begin(), added.end());
              matched.motifPath = entry.motifPath;
              matched.motifPath.push_back(motif.id);

              updateState(next,
                          trim(extended, frontier, options.globalInjectivity),
                          std::move(matched));
            });
      }
    }

    // Beam search, and the pathological-input guard, both prune here.
    const unsigned limit =
        options.beamWidth != 0
            ? options.beamWidth
            : (options.maxStatesPerStep != 0 ? options.maxStatesPerStep : 0);
    if (limit != 0 && next.size() > limit) {
      llvm::SmallVector<std::pair<StateKey, StateEntry>, 8> sorted(
          next.begin(), next.end());
      llvm::stable_sort(sorted, [](const auto &a, const auto &b) {
        return isBetter(a.second, b.second);
      });
      sorted.truncate(limit);

      StepTable pruned;
      for (auto &kv : sorted)
        pruned.insert(kv);
      next = std::move(pruned);
      result.beamTruncated = true;
    }
  }

  // argmax over the final table.
  const StepTable &finalTable = tables[m];
  if (finalTable.empty())
    return result;

  unsigned bestIndex = 0;
  bool haveBest = false;
  for (auto [i, item] : llvm::enumerate(finalTable)) {
    if (!haveBest || isBetter(item.second, finalTable.begin()[bestIndex].second)) {
      bestIndex = static_cast<unsigned>(i);
      haveBest = true;
    }
  }

  result.score = finalTable.begin()[bestIndex].second.score;

  // Walk the back-pointers, recording each step's motif and the bindings it
  // introduced. Steps are visited last to first, so the reaction at step t is
  // prep.order[t].
  unsigned index = bestIndex;
  for (unsigned t = m; t-- > 0;) {
    const StateEntry &entry = tables[t + 1].begin()[index].second;
    if (entry.motif)
      result.matches[prep.order[t]] = *entry.motif;
    for (const Binding &b : entry.introduced)
      result.assignment[b.species] = b.concrete;
    index = entry.previous;
  }

  for (unsigned r = 0; r < query.numReactions; ++r)
    if (!result.matches[r])
      result.unmatched.push_back(r);

  // Report any chemical bound to two abstract species. Impossible when global
  // injectivity is enforced; possible otherwise, because trim() forgets a
  // binding once its species leaves the frontier.
  llvm::SmallVector<std::pair<SpeciesId, AbstractSpecies>, 8> used;
  for (auto [s, bound] : llvm::enumerate(result.assignment))
    if (bound)
      used.push_back({*bound, static_cast<AbstractSpecies>(s)});
  llvm::sort(used);
  for (size_t i = 1; i < used.size(); ++i)
    if (used[i].first == used[i - 1].first &&
        !llvm::is_contained(result.injectivityViolations, used[i].first))
      result.injectivityViolations.push_back(used[i].first);

  return result;
}
