//===- MotifDatabaseFixture.cpp - Compiled-in motif database ----*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Loads the TableGen-generated fixture, so lit tests run without any external
// data file.
//
//===----------------------------------------------------------------------===//

#include "Standalone/CRN2DB/MotifDatabaseLoader.h"

#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/SmallVector.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/Support/Error.h"

#include <algorithm>
#include <cstdint>
#include <initializer_list>
#include <limits>

#define GET_MOTIF_DB
#include "Standalone/CRN2DB/MotifsFixture.inc"

using namespace chemcomp;

namespace {

llvm::Error makeError(const llvm::Twine &msg) {
  return llvm::createStringError(llvm::inconvertibleErrorCode(), msg);
}

/// Add one generated row to `db`.
void addFixtureMotif(const FixtureMotif &row, MotifDatabase &db) {
  llvm::SmallVector<llvm::StringRef, 8> species(row.species.begin(),
                                                row.species.end());
  llvm::SmallVector<int32_t, 8> stoich(row.stoich.begin(), row.stoich.end());

  llvm::SmallVector<double, 2> rateRange;
  if (row.hasRate) {
    rateRange.push_back(row.rateLo);
    rateRange.push_back(row.rateHi);
  }

  // The generator already rejected length mismatches and repeated species, so
  // this cannot fail; assert rather than propagate an impossible error.
  bool added = db.addMotif(row.id, species, stoich, row.deltaG, row.deltaGErr,
                           rateRange, row.source);
  assert(added && "generated fixture row rejected by addMotif");
  (void)added;
}

const FixtureSlice *findSlice(llvm::StringRef name) {
  for (const FixtureSlice &slice : kFixtureSlices)
    if (name == slice.name)
      return &slice;
  return nullptr;
}

} // namespace

llvm::SmallVector<llvm::StringRef, 8> chemcomp::getFixtureSliceNames() {
  llvm::SmallVector<llvm::StringRef, 8> names;
  for (const FixtureSlice &slice : kFixtureSlices)
    names.push_back(slice.name);
  // Already sorted by the generator; sorting again keeps that an invariant of
  // this function rather than of the generated file.
  llvm::sort(names);
  return names;
}

llvm::Error chemcomp::loadFixtureDatabase(llvm::StringRef sliceName,
                                          const SliceSpec &spec,
                                          MotifDatabase &db) {
  if (sliceName.empty()) {
    for (const FixtureMotif &row : kFixtureMotifs)
      addFixtureMotif(row, db);
  } else {
    const FixtureSlice *slice = findSlice(sliceName);
    if (!slice) {
      // An unknown slice would otherwise present as an empty database, which
      // looks like "nothing matched" rather than "you typed it wrong".
      std::string available;
      for (auto [i, name] : llvm::enumerate(getFixtureSliceNames())) {
        if (i)
          available += ", ";
        available += name.str();
      }
      return makeError("unknown fixture slice '" + sliceName +
                       "'; available slices: " + available);
    }

    // Membership is stored as ids; resolve against the motif table. Both are
    // sorted by id, so this stays linear.
    for (const char *id : slice->motifIds) {
      const FixtureMotif *row = llvm::find_if(
          kFixtureMotifs,
          [&](const FixtureMotif &m) { return llvm::StringRef(m.id) == id; });
      assert(row != std::end(kFixtureMotifs) &&
             "slice references an unknown motif id");
      if (row != std::end(kFixtureMotifs))
        addFixtureMotif(*row, db);
    }
  }

  db.applySlice(spec);
  db.buildIndex();
  return llvm::Error::success();
}
