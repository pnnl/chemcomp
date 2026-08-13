//===- MotifDatabase.h - Reaction motif database ----------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// In-memory representation of a chemical reaction motif database, as consumed
// by the CRN2DB matching pass (see papers/cf26-crn2db).
//
// The structure is deliberately loader-agnostic: it is populated identically
// from the TableGen-generated compiled-in fixture and from a runtime JSON
// database, so the matching algorithm never learns which one it got.
//
// Determinism is a hard requirement here. Match results are checked by
// FileCheck, so anything the algorithm iterates over must have a stable order
// independent of hash seeds or pointer values. Concretely: motifs are sorted by
// id, every candidate/species list is a sorted vector, and no observable
// iteration ever runs over a DenseMap/StringMap directly.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASE_H
#define CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASE_H

#include "llvm/ADT/ArrayRef.h"
#include "llvm/ADT/DenseMap.h"
#include "llvm/ADT/DenseMapInfo.h"
#include "llvm/ADT/Hashing.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringMap.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/Support/Allocator.h"
#include "llvm/Support/StringSaver.h"
#include "llvm/Support/raw_ostream.h"

#include <cmath>
#include <cstdint>
#include <limits>

namespace chemcomp {

/// Dense index into MotifDatabase's species table.
using SpeciesId = unsigned;
/// Dense index into MotifDatabase's motif table.
using MotifIdx = unsigned;

/// Sentinel used by ModelSEED for "no thermodynamic data". 43.1% of rows in
/// reactions_reduced.csv carry it, so it must never be mistaken for a value.
constexpr double kUnknownDeltaGSentinel = 1e7;

/// True when a deltaG/deltaGerr reading carries no information.
inline bool isUnknownDeltaG(double v) {
  return std::isnan(v) || v >= kUnknownDeltaGSentinel;
}

//===----------------------------------------------------------------------===//
// SignatureKey
//===----------------------------------------------------------------------===//

/// The sorted multiset of non-zero stoichiometric coefficients of a reaction,
/// e.g. `(-2,-1,2)`. Paper section 4.4 keys the candidate index on this.
///
/// The payload is an ArrayRef into arena-allocated storage owned by the
/// MotifDatabase, so a SignatureKey is cheap to copy but only valid for the
/// lifetime of the database that produced it.
struct SignatureKey {
  llvm::ArrayRef<int32_t> coeffs;

  SignatureKey() = default;
  explicit SignatureKey(llvm::ArrayRef<int32_t> c) : coeffs(c) {}

  bool operator==(const SignatureKey &o) const { return coeffs == o.coeffs; }
  bool operator!=(const SignatureKey &o) const { return !(*this == o); }

  /// Lexicographic order, for deterministic iteration over signatures.
  bool operator<(const SignatureKey &o) const {
    return std::lexicographical_compare(coeffs.begin(), coeffs.end(),
                                        o.coeffs.begin(), o.coeffs.end());
  }

  /// True when every element of `this` appears in `o` with at least the same
  /// multiplicity. Both must be sorted ascending.
  ///
  /// This implements the paper's coverage rule `phi(supp(r)) subseteq D_j`:
  /// a motif carrying extra species legitimately matches a smaller query, but
  /// has a longer signature, so exact-key lookup alone would miss it.
  bool isSubMultisetOf(const SignatureKey &o) const;

  void print(llvm::raw_ostream &os) const;
};

//===----------------------------------------------------------------------===//
// Motif
//===----------------------------------------------------------------------===//

/// One concrete reaction from the database.
///
/// `species` is sorted ascending and duplicate-free (the "set" of the paper's
/// `D_j`); `stoich` is positionally parallel to it, together forming the map
/// `nu_j : D_j -> Z`. Negative entries are reactants, positive are products,
/// matching ChemComp's convention throughout.
struct Motif {
  llvm::StringRef id;
  llvm::SmallVector<SpeciesId, 8> species;
  llvm::SmallVector<int32_t, 8> stoich;

  /// NaN when unknown. Parsed but unused by the matcher today: the paper sets
  /// the thermodynamic weights w_g = w_u = 0.
  double deltaG = std::numeric_limits<double>::quiet_NaN();
  double deltaGErr = std::numeric_limits<double>::quiet_NaN();

  /// Rate bounds. The database schema has no scalar rate constant, so the
  /// nominal kappa_j is derived from these; see `getNominalRate()`.
  double rateLo = 0.0;
  double rateHi = 0.0;
  bool hasRate = false;

  llvm::StringRef source;

  /// Sorted multiset of non-zero coefficients. Filled by buildIndex().
  SignatureKey signature;

  /// Stoichiometric coefficient of `s`, or 0 when `s` is absent.
  /// O(log n) via binary search over the sorted `species`.
  int32_t stoichFor(SpeciesId s) const;

  bool hasSpecies(SpeciesId s) const;

  /// Nominal rate kappa_j, as the geometric mean of the rate range.
  ///
  /// The geometric mean is the right center because S_r compares rates in log
  /// space. Returns NaN when no usable range is present, which the scorer
  /// treats as "unknown" and scores as 1.0 -- the same convention the paper
  /// uses for an unknown query rate k_i.
  double getNominalRate() const;
};

//===----------------------------------------------------------------------===//
// SliceSpec
//===----------------------------------------------------------------------===//

/// A declarative subset filter applied at load time.
///
/// Slicing is not merely a test-runtime convenience. Measured over the full
/// ModelSEED database (43,455 integer-coefficient motifs), only 913 distinct
/// signatures exist, and the single bucket `(-1,-1,1,1)` holds 11,421 motifs;
/// 90.2% of motifs sit in buckets larger than 100. Signature indexing alone
/// therefore does not bound the candidate set the way the paper assumes, and
/// filtering is what keeps the DP tractable.
struct SliceSpec {
  /// Keep a motif when it contains at least one of these species (by name).
  /// Empty means no species filter.
  llvm::SmallVector<std::string, 4> requireAnySpecies;

  /// Keep a motif whose `source` starts with this prefix, e.g. "ModelSEED:".
  /// Empty means no source filter.
  std::string sourcePrefix;

  /// Cap on motifs retained per signature bucket; 0 means unlimited.
  /// Truncation keeps the lowest motif ids, so it is deterministic.
  unsigned maxPerSignature = 0;

  /// Cap on total motifs retained; 0 means unlimited.
  unsigned maxMotifs = 0;

  bool isTrivial() const {
    return requireAnySpecies.empty() && sourcePrefix.empty() &&
           maxPerSignature == 0 && maxMotifs == 0;
  }
};

//===----------------------------------------------------------------------===//
// MotifDatabase
//===----------------------------------------------------------------------===//

/// A queryable collection of reaction motifs.
///
/// Typical use: call addMotif() repeatedly from a loader, then buildIndex()
/// once, then query. Queries before buildIndex() are not valid.
class MotifDatabase {
public:
  MotifDatabase() = default;

  // Non-copyable: StringRefs and SignatureKeys point into `alloc`.
  MotifDatabase(const MotifDatabase &) = delete;
  MotifDatabase &operator=(const MotifDatabase &) = delete;

  //===--------------------------------------------------------------------===//
  // Construction
  //===--------------------------------------------------------------------===//

  /// Intern a concrete species name, assigning a fresh id on first sight.
  SpeciesId internSpecies(llvm::StringRef name);

  /// Add one motif. `speciesNames` and `stoich` must be the same length and
  /// positionally parallel; order need not be sorted and is normalized here.
  /// Zero coefficients are dropped. Returns failure (false) on a length
  /// mismatch or a duplicate species within the motif.
  bool addMotif(llvm::StringRef id, llvm::ArrayRef<llvm::StringRef> speciesNames,
                llvm::ArrayRef<int32_t> stoich, double deltaG, double deltaGErr,
                llvm::ArrayRef<double> rateRange, llvm::StringRef source);

  /// Sort motifs by id, compute signatures, and build the lookup indices.
  /// Must be called once after all addMotif() calls and before any query.
  void buildIndex();

  /// Drop motifs failing `spec`. Call before buildIndex().
  /// Returns the number of motifs removed.
  unsigned applySlice(const SliceSpec &spec);

  //===--------------------------------------------------------------------===//
  // Queries
  //===--------------------------------------------------------------------===//

  size_t getNumMotifs() const { return motifs.size(); }
  size_t getNumSpecies() const { return idToName.size(); }

  const Motif &getMotif(MotifIdx i) const { return motifs[i]; }
  llvm::ArrayRef<Motif> getMotifs() const { return motifs; }

  llvm::StringRef getSpeciesName(SpeciesId s) const { return idToName[s]; }

  /// Look up a species by name; returns false when absent.
  bool lookupSpecies(llvm::StringRef name, SpeciesId &out) const;

  /// Motifs whose signature is exactly `sig`, in ascending motif-id order.
  llvm::ArrayRef<MotifIdx> lookupExact(const SignatureKey &sig) const;

  /// Motifs whose signature contains `sig` as a sub-multiset, in ascending
  /// motif-id order. This is the coverage-correct lookup; `lookupExact` is the
  /// paper-literal one. Results are memoized per distinct query signature.
  llvm::ArrayRef<MotifIdx> lookupSupersets(const SignatureKey &sig) const;

  /// Motifs containing `s`, in ascending motif-id order. Drives species-based
  /// slicing and seeds AC-3 domain construction.
  llvm::ArrayRef<MotifIdx> lookupBySpecies(SpeciesId s) const;

  /// All distinct signatures present, sorted. Small in practice: 913 over the
  /// full ModelSEED database.
  llvm::ArrayRef<SignatureKey> getSignatures() const { return signatureKeys; }

  /// Build the signature of an arbitrary coefficient row, interning the sorted
  /// key into this database's arena. Zeros are dropped.
  SignatureKey makeSignature(llvm::ArrayRef<int32_t> coeffs);

  //===--------------------------------------------------------------------===//
  // Debugging
  //===--------------------------------------------------------------------===//

  /// Stable textual dump used by lit tests to compare loaders. The format is
  /// deliberately compact and fully ordered.
  void dump(llvm::raw_ostream &os) const;

private:
  llvm::BumpPtrAllocator alloc;
  /// Interns names so the 140 duplicated compound names in the real database
  /// cost one string each, and every StringRef stays valid for our lifetime.
  llvm::UniqueStringSaver saver{alloc};

  llvm::StringMap<SpeciesId> nameToId;
  llvm::SmallVector<llvm::StringRef, 0> idToName;

  /// Sorted by id after buildIndex(), so MotifIdx order == motif-id order and
  /// every candidate list built from it is deterministic.
  llvm::SmallVector<Motif, 0> motifs;

  llvm::DenseMap<SignatureKey, llvm::SmallVector<MotifIdx, 8>> index;
  llvm::SmallVector<SignatureKey, 0> signatureKeys;
  llvm::DenseMap<SpeciesId, llvm::SmallVector<MotifIdx, 8>> bySpecies;

  /// Memo for lookupSupersets(). Mutable: it is a pure cache.
  mutable llvm::DenseMap<SignatureKey, llvm::SmallVector<MotifIdx, 8>>
      supersetCache;

  bool indexBuilt = false;

  /// Copy `coeffs` into the arena and wrap it as a key.
  SignatureKey internSignature(llvm::ArrayRef<int32_t> coeffs);
};

} // namespace chemcomp

//===----------------------------------------------------------------------===//
// DenseMapInfo for SignatureKey
//===----------------------------------------------------------------------===//

namespace llvm {
template <> struct DenseMapInfo<chemcomp::SignatureKey> {
  static inline chemcomp::SignatureKey getEmptyKey() {
    return chemcomp::SignatureKey(ArrayRef<int32_t>(
        reinterpret_cast<const int32_t *>(~uintptr_t(0)), size_t(0)));
  }
  static inline chemcomp::SignatureKey getTombstoneKey() {
    return chemcomp::SignatureKey(ArrayRef<int32_t>(
        reinterpret_cast<const int32_t *>(~uintptr_t(1)), size_t(0)));
  }
  static unsigned getHashValue(const chemcomp::SignatureKey &k) {
    return static_cast<unsigned>(
        hash_combine_range(k.coeffs.begin(), k.coeffs.end()));
  }
  static bool isEqual(const chemcomp::SignatureKey &a,
                      const chemcomp::SignatureKey &b) {
    // Empty/tombstone carry sentinel pointers that must not be dereferenced,
    // so compare identity before touching the payload.
    if (a.coeffs.data() == b.coeffs.data() &&
        a.coeffs.size() == b.coeffs.size())
      return true;
    if (a.coeffs.data() == getEmptyKey().coeffs.data() ||
        a.coeffs.data() == getTombstoneKey().coeffs.data() ||
        b.coeffs.data() == getEmptyKey().coeffs.data() ||
        b.coeffs.data() == getTombstoneKey().coeffs.data())
      return false;
    return a == b;
  }
};
} // namespace llvm

#endif // CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASE_H
