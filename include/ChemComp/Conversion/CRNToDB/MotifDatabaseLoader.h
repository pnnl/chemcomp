//===- MotifDatabaseLoader.h - Motif database loaders -----------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// The two ways to populate a MotifDatabase.
//
// Both apply the same SliceSpec and call buildIndex(), so the resulting
// database is indistinguishable to the matcher regardless of origin -- a
// property asserted directly by test/Conversion/CRNToDB/db-loader-equivalence.mlir.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASELOADER_H
#define CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASELOADER_H

#include "ChemComp/Conversion/CRNToDB/MotifDatabase.h"

#include "llvm/ADT/StringRef.h"
#include "llvm/Support/Error.h"

namespace mlir::cc {

/// Populate `db` from the compiled-in TableGen fixture.
///
/// `sliceName` selects a slice declared in MotifsFixture.td; an empty name
/// means every fixture motif. Returns an error naming the available slices
/// when `sliceName` is unknown, since a typo would otherwise look like an
/// empty database.
llvm::Error loadFixtureDatabase(llvm::StringRef sliceName,
                                const SliceSpec &spec, MotifDatabase &db);

/// Names of the slices available to loadFixtureDatabase(), in sorted order.
llvm::SmallVector<llvm::StringRef, 8> getFixtureSliceNames();

/// Populate `db` from a JSON motif database.
///
/// The expected shape is a top-level array of objects:
///
///   [{"motif_id": "rxn00006",
///     "species": ["H2O2", "H2O", "O2"],
///     "stoich": [-2, 2, 1],
///     "deltaG": -46.06,          // optional; null or absent means unknown
///     "deltaGerr": 1.64,         // optional
///     "rate_range": [0.01, 1.0], // optional
///     "source": "KEGG:R00009|ModelSEED:rxn00006"}]
///
/// `species` and `stoich` must be the same length. Errors carry the offending
/// record's index and id so a bad entry in a 43k-record file is findable.
llvm::Error loadJSONDatabase(llvm::StringRef path, const SliceSpec &spec,
                             MotifDatabase &db);

/// Parse a JSON motif database from memory. Shares all validation with
/// loadJSONDatabase(); `context` labels diagnostics.
llvm::Error loadJSONDatabaseFromBuffer(llvm::StringRef json,
                                       llvm::StringRef context,
                                       const SliceSpec &spec,
                                       MotifDatabase &db);

} // namespace mlir::cc

#endif // CHEMCOMP_CONVERSION_CRNTODB_MOTIFDATABASELOADER_H
