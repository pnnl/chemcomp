//===- Preprocessing.cpp - CRN matching preprocessing -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/CRNToDB/CRNMatcher.h"

#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallSet.h"

#include <algorithm>

using namespace mlir::cc;

//===----------------------------------------------------------------------===//
// QueryCRN
//===----------------------------------------------------------------------===//

llvm::SmallVector<AbstractSpecies, 8>
QueryCRN::support(AbstractReaction r) const {
  llvm::SmallVector<AbstractSpecies, 8> out;
  for (auto [s, v] : llvm::enumerate(rows[r]))
    if (v != 0)
      out.push_back(static_cast<AbstractSpecies>(s));
  return out;
}

llvm::SmallVector<int32_t, 8>
QueryCRN::signatureCoeffs(AbstractReaction r) const {
  llvm::SmallVector<int32_t, 8> out;
  for (int32_t v : rows[r])
    if (v != 0)
      out.push_back(v);
  llvm::sort(out);
  return out;
}

//===----------------------------------------------------------------------===//
// Min-fill elimination ordering
//===----------------------------------------------------------------------===//

namespace {

/// Dense adjacency over abstract species.
struct InteractionGraph {
  unsigned n = 0;
  llvm::SmallVector<llvm::SmallBitVector, 8> adjacency;

  explicit InteractionGraph(const QueryCRN &query) : n(query.numSpecies) {
    adjacency.assign(n, llvm::SmallBitVector(n, false));
    // Two species interact when some reaction touches both: the reaction is a
    // factor coupling them, in the paper's factor-graph reading.
    for (unsigned r = 0; r < query.numReactions; ++r) {
      llvm::SmallVector<AbstractSpecies, 8> supp = query.support(r);
      for (size_t i = 0; i < supp.size(); ++i)
        for (size_t j = i + 1; j < supp.size(); ++j) {
          adjacency[supp[i]].set(supp[j]);
          adjacency[supp[j]].set(supp[i]);
        }
    }
  }

  /// Fill edges added by eliminating `v`: pairs of live neighbours not already
  /// adjacent.
  unsigned fillCount(AbstractSpecies v, const llvm::SmallBitVector &live) const {
    llvm::SmallVector<AbstractSpecies, 8> neighbours;
    for (unsigned u = 0; u < n; ++u)
      if (live[u] && u != v && adjacency[v][u])
        neighbours.push_back(u);

    unsigned fill = 0;
    for (size_t i = 0; i < neighbours.size(); ++i)
      for (size_t j = i + 1; j < neighbours.size(); ++j)
        if (!adjacency[neighbours[i]][neighbours[j]])
          ++fill;
    return fill;
  }

  unsigned degree(AbstractSpecies v, const llvm::SmallBitVector &live) const {
    unsigned d = 0;
    for (unsigned u = 0; u < n; ++u)
      if (live[u] && u != v && adjacency[v][u])
        ++d;
    return d;
  }

  /// Connect `v`'s neighbours to each other, then drop `v`.
  void eliminate(AbstractSpecies v, llvm::SmallBitVector &live) {
    llvm::SmallVector<AbstractSpecies, 8> neighbours;
    for (unsigned u = 0; u < n; ++u)
      if (live[u] && u != v && adjacency[v][u])
        neighbours.push_back(u);

    for (size_t i = 0; i < neighbours.size(); ++i)
      for (size_t j = i + 1; j < neighbours.size(); ++j) {
        adjacency[neighbours[i]].set(neighbours[j]);
        adjacency[neighbours[j]].set(neighbours[i]);
      }
    live.reset(v);
  }
};

} // namespace

llvm::SmallVector<AbstractSpecies, 8>
mlir::cc::computeMinFillOrder(const QueryCRN &query) {
  InteractionGraph graph(query);
  llvm::SmallBitVector live(query.numSpecies, true);
  llvm::SmallVector<AbstractSpecies, 8> order;
  order.reserve(query.numSpecies);

  for (unsigned step = 0; step < query.numSpecies; ++step) {
    AbstractSpecies best = 0;
    unsigned bestFill = std::numeric_limits<unsigned>::max();
    unsigned bestDegree = std::numeric_limits<unsigned>::max();
    bool found = false;

    for (unsigned v = 0; v < query.numSpecies; ++v) {
      if (!live[v])
        continue;
      unsigned fill = graph.fillCount(v, live);
      unsigned degree = graph.degree(v, live);
      // Fewest fill edges wins; then lowest degree; then lowest index, so the
      // ordering never depends on iteration order.
      if (!found || fill < bestFill ||
          (fill == bestFill && degree < bestDegree)) {
        best = v;
        bestFill = fill;
        bestDegree = degree;
        found = true;
      }
    }

    if (!found)
      break;
    order.push_back(best);
    graph.eliminate(best, live);
  }

  return order;
}

llvm::SmallVector<AbstractReaction, 8>
mlir::cc::deriveReactionOrder(const QueryCRN &query,
                              llvm::ArrayRef<AbstractSpecies> speciesOrder) {
  // position[s] = when s is eliminated.
  llvm::SmallVector<unsigned, 8> position(query.numSpecies, 0);
  for (auto [i, s] : llvm::enumerate(speciesOrder))
    position[s] = static_cast<unsigned>(i);

  struct Key {
    unsigned last;                                 // latest position in supp
    llvm::SmallVector<unsigned, 8> sortedPositions; // tie-break
    AbstractReaction reaction;                     // final tie-break
  };

  llvm::SmallVector<Key, 8> keys;
  keys.reserve(query.numReactions);
  for (unsigned r = 0; r < query.numReactions; ++r) {
    Key key;
    key.reaction = r;
    key.last = 0;
    for (AbstractSpecies s : query.support(r)) {
      key.sortedPositions.push_back(position[s]);
      key.last = std::max(key.last, position[s]);
    }
    llvm::sort(key.sortedPositions);
    keys.push_back(std::move(key));
  }

  llvm::sort(keys, [](const Key &a, const Key &b) {
    if (a.last != b.last)
      return a.last < b.last;
    if (a.sortedPositions != b.sortedPositions)
      return std::lexicographical_compare(
          a.sortedPositions.begin(), a.sortedPositions.end(),
          b.sortedPositions.begin(), b.sortedPositions.end());
    return a.reaction < b.reaction;
  });

  llvm::SmallVector<AbstractReaction, 8> order;
  order.reserve(keys.size());
  for (const Key &key : keys)
    order.push_back(key.reaction);
  return order;
}

llvm::SmallVector<llvm::SmallVector<AbstractSpecies, 8>, 8>
mlir::cc::computeFrontiers(const QueryCRN &query,
                           llvm::ArrayRef<AbstractReaction> order) {
  const unsigned m = static_cast<unsigned>(order.size());

  // suffix[t] = species appearing in reactions order[t..m-1].
  llvm::SmallVector<llvm::SmallBitVector, 8> suffix(m + 1,
                                                    llvm::SmallBitVector(
                                                        query.numSpecies,
                                                        false));
  for (unsigned t = m; t-- > 0;) {
    suffix[t] = suffix[t + 1];
    for (AbstractSpecies s : query.support(order[t]))
      suffix[t].set(s);
  }

  llvm::SmallVector<llvm::SmallVector<AbstractSpecies, 8>, 8> frontiers;
  frontiers.reserve(m + 1);

  llvm::SmallBitVector seen(query.numSpecies, false);
  // F_0 is empty: nothing has been processed yet.
  frontiers.push_back({});

  for (unsigned t = 0; t < m; ++t) {
    for (AbstractSpecies s : query.support(order[t]))
      seen.set(s);

    // F_t = (species seen so far) intersect (species still to come).
    llvm::SmallVector<AbstractSpecies, 8> live;
    for (unsigned s = 0; s < query.numSpecies; ++s)
      if (seen[s] && suffix[t + 1][s])
        live.push_back(s);
    frontiers.push_back(std::move(live));
  }

  return frontiers;
}

//===----------------------------------------------------------------------===//
// Preprocessing
//===----------------------------------------------------------------------===//

unsigned Preprocessing::maxFrontierWidth() const {
  unsigned width = 0;
  for (const auto &f : frontier)
    width = std::max(width, static_cast<unsigned>(f.size()));
  return width;
}

namespace {

/// Motifs that could serve `reaction`, before AC-3.
llvm::SmallVector<MotifIdx, 8>
collectCandidates(const QueryCRN &query, AbstractReaction reaction,
                  const MotifDatabase &db, const PreprocessOptions &options,
                  unsigned &untruncatedCount) {
  llvm::SmallVector<int32_t, 8> coeffs = query.signatureCoeffs(reaction);
  untruncatedCount = 0;
  if (coeffs.empty())
    return {};

  SignatureKey key(coeffs);
  llvm::ArrayRef<MotifIdx> found = options.strictSignature
                                       ? db.lookupExact(key)
                                       : db.lookupSupersets(key);
  untruncatedCount = static_cast<unsigned>(found.size());

  llvm::SmallVector<MotifIdx, 8> candidates(found.begin(), found.end());
  if (options.maxCandidates != 0 && candidates.size() > options.maxCandidates)
    candidates.truncate(options.maxCandidates);
  return candidates;
}

/// Concrete species that could bind `s`, given the reactions it participates in.
///
/// A concrete species is admissible only if, for every reaction touching `s`,
/// some candidate motif of that reaction offers it with the exact coefficient
/// the query demands. Intersecting across reactions is already a strong filter;
/// AC-3 then removes what survives here but cannot be part of any consistent
/// whole-reaction assignment.
llvm::SmallVector<SpeciesId, 8>
initialDomain(const QueryCRN &query, AbstractSpecies s, const MotifDatabase &db,
              llvm::ArrayRef<llvm::SmallVector<MotifIdx, 8>> candidates) {
  llvm::SmallVector<SpeciesId, 8> domain;
  bool first = true;

  for (unsigned r = 0; r < query.numReactions; ++r) {
    int32_t want = query.coeff(r, s);
    if (want == 0)
      continue;

    llvm::SmallVector<SpeciesId, 8> allowed;
    for (MotifIdx idx : candidates[r]) {
      const Motif &motif = db.getMotif(idx);
      for (auto [i, concrete] : llvm::enumerate(motif.species))
        if (motif.stoich[i] == want)
          allowed.push_back(concrete);
    }
    llvm::sort(allowed);
    allowed.erase(std::unique(allowed.begin(), allowed.end()), allowed.end());

    if (first) {
      domain = std::move(allowed);
      first = false;
      continue;
    }

    llvm::SmallVector<SpeciesId, 8> intersection;
    std::set_intersection(domain.begin(), domain.end(), allowed.begin(),
                          allowed.end(), std::back_inserter(intersection));
    domain = std::move(intersection);
    if (domain.empty())
      break;
  }

  return domain;
}

/// True when some candidate motif of `r` can bind `s` to `concrete` while also
/// offering a distinct concrete species for `other`.
///
/// This is the arc-consistency predicate: it asks whether the value `concrete`
/// has any support at all along the edge (s, other), which is what AC-3 prunes
/// on.
bool hasSupport(const QueryCRN &query, AbstractReaction r, AbstractSpecies s,
                SpeciesId concrete, AbstractSpecies other,
                llvm::ArrayRef<SpeciesId> otherDomain, const MotifDatabase &db,
                llvm::ArrayRef<MotifIdx> candidates) {
  int32_t wantS = query.coeff(r, s);
  int32_t wantOther = query.coeff(r, other);

  for (MotifIdx idx : candidates) {
    const Motif &motif = db.getMotif(idx);
    if (motif.stoichFor(concrete) != wantS)
      continue;

    for (auto [i, candidate] : llvm::enumerate(motif.species)) {
      if (motif.stoich[i] != wantOther)
        continue;
      // Injectivity: two abstract species cannot share one chemical.
      if (candidate == concrete)
        continue;
      if (std::binary_search(otherDomain.begin(), otherDomain.end(), candidate))
        return true;
    }
  }
  return false;
}

/// AC-3 over the species-sharing constraint graph.
void runArcConsistency(const QueryCRN &query, const MotifDatabase &db,
                       const PreprocessOptions &options, Preprocessing &result) {
  // Arcs connect species co-occurring in a reaction; that reaction carries the
  // constraint.
  struct Arc {
    AbstractSpecies from;
    AbstractSpecies to;
    AbstractReaction reaction;
  };

  llvm::SmallVector<Arc, 16> worklist;
  for (unsigned r = 0; r < query.numReactions; ++r) {
    llvm::SmallVector<AbstractSpecies, 8> supp = query.support(r);
    for (AbstractSpecies a : supp)
      for (AbstractSpecies b : supp)
        if (a != b)
          worklist.push_back({a, b, r});
  }

  unsigned steps = 0;
  while (!worklist.empty()) {
    if (++steps > options.arcConsistencyBudget) {
      // The worst case is O(|DB| * |Sigma|^2); rather than risk a hang, stop
      // and report that domains are larger than a full run would leave them.
      result.arcConsistencyIncomplete = true;
      return;
    }

    Arc arc = worklist.pop_back_val();
    llvm::SmallVector<SpeciesId, 8> &domain = result.domains[arc.from];
    if (domain.empty())
      continue;

    llvm::SmallVector<SpeciesId, 8> kept;
    for (SpeciesId concrete : domain)
      if (hasSupport(query, arc.reaction, arc.from, concrete, arc.to,
                     result.domains[arc.to], db, result.candidates[arc.reaction]))
        kept.push_back(concrete);

    if (kept.size() == domain.size())
      continue;

    domain = std::move(kept);

    // Shrinking a domain can invalidate support for its neighbours, so requeue
    // the arcs pointing back at it.
    for (unsigned r = 0; r < query.numReactions; ++r) {
      if (query.coeff(r, arc.from) == 0)
        continue;
      for (AbstractSpecies other : query.support(r))
        if (other != arc.from)
          worklist.push_back({other, arc.from, r});
    }
  }
}

} // namespace

Preprocessing mlir::cc::preprocess(const QueryCRN &query,
                                   const MotifDatabase &db,
                                   const PreprocessOptions &options) {
  Preprocessing result;

  result.order = deriveReactionOrder(query, computeMinFillOrder(query));
  result.frontier = computeFrontiers(query, result.order);

  result.candidates.resize(query.numReactions);
  for (unsigned r = 0; r < query.numReactions; ++r) {
    unsigned untruncated = 0;
    result.candidates[r] = collectCandidates(query, r, db, options, untruncated);
    if (untruncated > result.candidates[r].size())
      result.truncated.push_back({r, untruncated});
  }

  result.domains.resize(query.numSpecies);
  for (unsigned s = 0; s < query.numSpecies; ++s)
    result.domains[s] = initialDomain(query, s, db, result.candidates);

  runArcConsistency(query, db, options, result);
  return result;
}

//===----------------------------------------------------------------------===//
// Debug printing
//===----------------------------------------------------------------------===//

void Preprocessing::print(llvm::raw_ostream &os, const QueryCRN &query,
                          const MotifDatabase &db) const {
  os << "order:";
  for (AbstractReaction r : order)
    os << ' ' << r;
  os << "\n";

  os << "max-frontier-width: " << maxFrontierWidth() << "\n";

  for (auto [t, live] : llvm::enumerate(frontier)) {
    os << "frontier[" << t << "]:";
    for (AbstractSpecies s : live)
      os << ' ' << s;
    os << "\n";
  }

  for (auto [s, domain] : llvm::enumerate(domains)) {
    os << "domain[" << s << "]: " << domain.size();
    for (SpeciesId c : domain)
      os << ' ' << db.getSpeciesName(c);
    os << "\n";
  }

  for (auto [r, list] : llvm::enumerate(candidates)) {
    os << "candidates[" << r << "]: " << list.size();
    for (MotifIdx idx : list)
      os << ' ' << db.getMotif(idx).id;
    os << "\n";
  }

  for (auto [reaction, before] : truncated)
    os << "truncated[" << reaction << "]: " << before << " -> "
       << candidates[reaction].size() << "\n";

  if (arcConsistencyIncomplete)
    os << "arc-consistency: incomplete\n";
}
