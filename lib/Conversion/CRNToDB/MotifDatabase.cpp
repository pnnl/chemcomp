//===- MotifDatabase.cpp - Reaction motif database --------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "ChemComp/Conversion/CRNToDB/MotifDatabase.h"

#include "llvm/ADT/STLExtras.h"
#include "llvm/Support/Format.h"

#include <algorithm>

using namespace chemcomp;

//===----------------------------------------------------------------------===//
// SignatureKey
//===----------------------------------------------------------------------===//

bool SignatureKey::isSubMultisetOf(const SignatureKey &o) const {
  // Both sides are sorted ascending, so a linear merge suffices.
  if (coeffs.size() > o.coeffs.size())
    return false;
  size_t i = 0, j = 0;
  while (i < coeffs.size() && j < o.coeffs.size()) {
    if (coeffs[i] == o.coeffs[j]) {
      ++i;
      ++j;
    } else if (o.coeffs[j] < coeffs[i]) {
      ++j;
    } else {
      // o[j] > this[i]: this[i] can no longer be matched.
      return false;
    }
  }
  return i == coeffs.size();
}

void SignatureKey::print(llvm::raw_ostream &os) const {
  os << '(';
  llvm::interleaveComma(coeffs, os);
  os << ')';
}

//===----------------------------------------------------------------------===//
// Motif
//===----------------------------------------------------------------------===//

int32_t Motif::stoichFor(SpeciesId s) const {
  auto it = std::lower_bound(species.begin(), species.end(), s);
  if (it == species.end() || *it != s)
    return 0;
  return stoich[std::distance(species.begin(), it)];
}

bool Motif::hasSpecies(SpeciesId s) const {
  return std::binary_search(species.begin(), species.end(), s);
}

double Motif::getNominalRate() const {
  if (!hasRate || rateLo <= 0.0 || rateHi <= 0.0)
    return std::numeric_limits<double>::quiet_NaN();
  return std::sqrt(rateLo * rateHi);
}

//===----------------------------------------------------------------------===//
// MotifDatabase - construction
//===----------------------------------------------------------------------===//

SpeciesId MotifDatabase::internSpecies(llvm::StringRef name) {
  auto it = nameToId.find(name);
  if (it != nameToId.end())
    return it->second;
  llvm::StringRef saved = saver.save(name);
  SpeciesId id = static_cast<SpeciesId>(idToName.size());
  nameToId[saved] = id;
  idToName.push_back(saved);
  return id;
}

bool MotifDatabase::lookupSpecies(llvm::StringRef name, SpeciesId &out) const {
  auto it = nameToId.find(name);
  if (it == nameToId.end())
    return false;
  out = it->second;
  return true;
}

SignatureKey MotifDatabase::internSignature(llvm::ArrayRef<int32_t> coeffs) {
  if (coeffs.empty())
    return SignatureKey(llvm::ArrayRef<int32_t>());
  auto *storage = alloc.Allocate<int32_t>(coeffs.size());
  std::copy(coeffs.begin(), coeffs.end(), storage);
  return SignatureKey(llvm::ArrayRef<int32_t>(storage, coeffs.size()));
}

SignatureKey MotifDatabase::makeSignature(llvm::ArrayRef<int32_t> coeffs) {
  llvm::SmallVector<int32_t, 8> nonZero;
  for (int32_t c : coeffs)
    if (c != 0)
      nonZero.push_back(c);
  llvm::sort(nonZero);
  return internSignature(nonZero);
}

bool MotifDatabase::addMotif(llvm::StringRef id,
                             llvm::ArrayRef<llvm::StringRef> speciesNames,
                             llvm::ArrayRef<int32_t> stoich, double deltaG,
                             double deltaGErr, llvm::ArrayRef<double> rateRange,
                             llvm::StringRef source) {
  if (speciesNames.size() != stoich.size())
    return false;

  // Pair up, dropping zero coefficients: a species with coefficient 0 does not
  // participate and would corrupt the signature.
  llvm::SmallVector<std::pair<SpeciesId, int32_t>, 8> entries;
  entries.reserve(speciesNames.size());
  for (auto [name, coeff] : llvm::zip_equal(speciesNames, stoich)) {
    if (coeff == 0)
      continue;
    entries.emplace_back(internSpecies(name), coeff);
  }

  llvm::sort(entries, [](const auto &a, const auto &b) {
    return a.first < b.first;
  });

  // Reject a species repeated within one motif: nu_j would be ambiguous.
  for (size_t i = 1; i < entries.size(); ++i)
    if (entries[i].first == entries[i - 1].first)
      return false;

  Motif m;
  m.id = saver.save(id);
  m.source = source.empty() ? llvm::StringRef() : saver.save(source);
  m.species.reserve(entries.size());
  m.stoich.reserve(entries.size());
  for (auto &e : entries) {
    m.species.push_back(e.first);
    m.stoich.push_back(e.second);
  }

  m.deltaG = isUnknownDeltaG(deltaG) ? std::numeric_limits<double>::quiet_NaN()
                                     : deltaG;
  m.deltaGErr = isUnknownDeltaG(deltaGErr)
                    ? std::numeric_limits<double>::quiet_NaN()
                    : deltaGErr;

  if (rateRange.size() == 2 && rateRange[0] > 0.0 && rateRange[1] > 0.0) {
    m.rateLo = std::min(rateRange[0], rateRange[1]);
    m.rateHi = std::max(rateRange[0], rateRange[1]);
    m.hasRate = true;
  }

  motifs.push_back(std::move(m));
  indexBuilt = false;
  return true;
}

unsigned MotifDatabase::applySlice(const SliceSpec &spec) {
  if (spec.isTrivial())
    return 0;

  const size_t before = motifs.size();

  // Resolve the species allowlist to ids once. A name absent from the database
  // simply matches nothing.
  llvm::SmallVector<SpeciesId, 4> required;
  for (const std::string &n : spec.requireAnySpecies) {
    SpeciesId id;
    if (lookupSpecies(n, id))
      required.push_back(id);
  }

  llvm::erase_if(motifs, [&](const Motif &m) {
    if (!spec.sourcePrefix.empty() && !m.source.starts_with(spec.sourcePrefix))
      return true;
    if (!spec.requireAnySpecies.empty()) {
      bool any = llvm::any_of(
          required, [&](SpeciesId s) { return m.hasSpecies(s); });
      if (!any)
        return true;
    }
    return false;
  });

  // Per-signature and global caps need a deterministic notion of "first", so
  // sort by id before truncating. buildIndex() re-sorts anyway; doing it here
  // makes the truncation reproducible regardless of insertion order.
  llvm::sort(motifs, [](const Motif &a, const Motif &b) { return a.id < b.id; });

  if (spec.maxPerSignature != 0) {
    // Count per signature as we scan in id order, dropping past the cap.
    llvm::DenseMap<SignatureKey, unsigned> seen;
    llvm::SmallVector<Motif, 0> kept;
    kept.reserve(motifs.size());
    for (Motif &m : motifs) {
      llvm::SmallVector<int32_t, 8> sorted(m.stoich.begin(), m.stoich.end());
      llvm::sort(sorted);
      SignatureKey sig = internSignature(sorted);
      unsigned &count = seen[sig];
      if (count < spec.maxPerSignature) {
        ++count;
        kept.push_back(std::move(m));
      }
    }
    motifs = std::move(kept);
  }

  if (spec.maxMotifs != 0 && motifs.size() > spec.maxMotifs)
    motifs.truncate(spec.maxMotifs);

  indexBuilt = false;
  return static_cast<unsigned>(before - motifs.size());
}

void MotifDatabase::buildIndex() {
  index.clear();
  signatureKeys.clear();
  bySpecies.clear();
  supersetCache.clear();

  // Sorting by id is what makes MotifIdx order == motif-id order, which in
  // turn makes every candidate list below deterministic without extra sorting.
  llvm::sort(motifs, [](const Motif &a, const Motif &b) { return a.id < b.id; });

  for (auto [i, m] : llvm::enumerate(motifs)) {
    llvm::SmallVector<int32_t, 8> sorted(m.stoich.begin(), m.stoich.end());
    llvm::sort(sorted);
    m.signature = internSignature(sorted);

    index[m.signature].push_back(static_cast<MotifIdx>(i));
    for (SpeciesId s : m.species)
      bySpecies[s].push_back(static_cast<MotifIdx>(i));
  }

  // Materialize the signature list in sorted order for deterministic scans.
  signatureKeys.reserve(index.size());
  for (const auto &kv : index)
    signatureKeys.push_back(kv.first);
  llvm::sort(signatureKeys);

  indexBuilt = true;
}

//===----------------------------------------------------------------------===//
// MotifDatabase - queries
//===----------------------------------------------------------------------===//

llvm::ArrayRef<MotifIdx>
MotifDatabase::lookupExact(const SignatureKey &sig) const {
  assert(indexBuilt && "lookupExact() before buildIndex()");
  auto it = index.find(sig);
  if (it == index.end())
    return {};
  return it->second;
}

llvm::ArrayRef<MotifIdx>
MotifDatabase::lookupSupersets(const SignatureKey &sig) const {
  assert(indexBuilt && "lookupSupersets() before buildIndex()");

  auto cached = supersetCache.find(sig);
  if (cached != supersetCache.end())
    return cached->second;

  // Scan the (small: ~913 in the full database) signature key set for keys
  // containing `sig`, then union their buckets. Because signatureKeys is
  // sorted and each bucket is in ascending MotifIdx order, the union is built
  // by merging and then sorting once -- cheap, and deterministic.
  llvm::SmallVector<MotifIdx, 8> result;
  for (const SignatureKey &key : signatureKeys) {
    if (!sig.isSubMultisetOf(key))
      continue;
    auto it = index.find(key);
    if (it != index.end())
      result.append(it->second.begin(), it->second.end());
  }
  llvm::sort(result);
  result.erase(std::unique(result.begin(), result.end()), result.end());

  auto inserted = supersetCache.try_emplace(sig, std::move(result));
  return inserted.first->second;
}

llvm::ArrayRef<MotifIdx> MotifDatabase::lookupBySpecies(SpeciesId s) const {
  assert(indexBuilt && "lookupBySpecies() before buildIndex()");
  auto it = bySpecies.find(s);
  if (it == bySpecies.end())
    return {};
  return it->second;
}

//===----------------------------------------------------------------------===//
// MotifDatabase - debugging
//===----------------------------------------------------------------------===//

void MotifDatabase::dump(llvm::raw_ostream &os) const {
  os << "motifs: " << motifs.size() << "\n";
  os << "species: " << idToName.size() << "\n";
  os << "signatures: " << signatureKeys.size() << "\n";

  for (const Motif &m : motifs) {
    os << "motif " << m.id << " sig=";
    m.signature.print(os);
    os << " {";
    for (auto [i, s] : llvm::enumerate(m.species)) {
      if (i)
        os << ", ";
      os << idToName[s] << ':' << m.stoich[i];
    }
    os << "}";

    if (std::isnan(m.deltaG))
      os << " deltaG=unknown";
    else
      os << llvm::format(" deltaG=%g", m.deltaG);

    if (m.hasRate)
      os << llvm::format(" rate=[%g,%g]", m.rateLo, m.rateHi);
    else
      os << " rate=unknown";

    if (!m.source.empty())
      os << " source=" << m.source;
    os << "\n";
  }
}
