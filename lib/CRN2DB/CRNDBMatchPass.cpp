//===- CRNDBMatchPass.cpp - CRN to database matching ------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Implements the `crn-db-match` pass: binds the abstract species and reactions
// of an `acr.crn_generic` operation to concrete database motifs.
//
// See papers/cf26-crn2db. The matching algorithm itself lands in a later
// change; this file currently handles database loading, query extraction, and
// result annotation, with a greedy stand-in for the dynamic program.
//
//===----------------------------------------------------------------------===//

#include "Standalone/CRN2DB/CRN2DBPasses.h"

#include "Standalone/CRN2DB/MotifDatabase.h"
#include "Standalone/CRN2DB/MotifDatabaseLoader.h"

#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/IR/BuiltinTypes.h"

#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/raw_ostream.h"

using namespace mlir;
using namespace chemcomp;

namespace mlir::crn2db {
#define GEN_PASS_DEF_CRNDBMATCH
#include "Standalone/CRN2DB/CRN2DBPasses.h.inc"

namespace {

/// The abstract CRN extracted from an `acr.crn_generic` operation: the paper's
/// query Q = (S_Q, R_Q, N_Q, k_Q), minus the rate constants for now.
struct QueryCRN {
  /// Row-major stoichiometry, reactions x species (ChemComp's convention).
  llvm::SmallVector<llvm::SmallVector<int32_t, 8>, 8> rows;
  int64_t numReactions = 0;
  int64_t numSpecies = 0;

  /// Column indices with a non-zero coefficient in row `r`: supp(r).
  llvm::SmallVector<unsigned, 8> support(unsigned r) const {
    llvm::SmallVector<unsigned, 8> out;
    for (auto [c, v] : llvm::enumerate(rows[r]))
      if (v != 0)
        out.push_back(static_cast<unsigned>(c));
    return out;
  }
};

/// Read the `stoich_matrix` attribute into a QueryCRN.
///
/// The operation verifier has already checked that the attribute exists, is
/// 2-D, and agrees with the operand counts, so failure here is not expected.
bool extractQuery(acr::CRNGenericOp op, QueryCRN &query) {
  auto stoichAttr = op->getAttrOfType<DenseIntElementsAttr>("stoich_matrix");
  if (!stoichAttr)
    return false;

  auto type = llvm::dyn_cast<ShapedType>(stoichAttr.getType());
  if (!type || type.getRank() != 2)
    return false;

  query.numReactions = type.getDimSize(0);
  query.numSpecies = type.getDimSize(1);

  auto values = stoichAttr.getValues<llvm::APInt>();
  auto it = values.begin();
  query.rows.reserve(query.numReactions);
  for (int64_t r = 0; r < query.numReactions; ++r) {
    llvm::SmallVector<int32_t, 8> row;
    row.reserve(query.numSpecies);
    for (int64_t c = 0; c < query.numSpecies; ++c, ++it)
      row.push_back(static_cast<int32_t>((*it).getSExtValue()));
    query.rows.push_back(std::move(row));
  }
  return true;
}

class CRNDBMatch : public impl::CRNDBMatchBase<CRNDBMatch> {
public:
  using impl::CRNDBMatchBase<CRNDBMatch>::CRNDBMatchBase;

  void runOnOperation() final {
    ModuleOp module = getOperation();

    // Build the slice from the pass options once; it applies to whichever
    // loader runs.
    SliceSpec spec;
    for (const std::string &s : filterSpecies)
      spec.requireAnySpecies.push_back(s);
    spec.sourcePrefix = filterSource;
    spec.maxPerSignature = maxCandidates;

    MotifDatabase database;
    if (llvm::Error error = loadDatabase(spec, database)) {
      module.emitError() << "crn-db-match: "
                         << llvm::toString(std::move(error));
      return signalPassFailure();
    }

    if (dumpDatabase) {
      database.dump(llvm::outs());
      return;
    }

    module.walk([&](acr::CRNGenericOp op) { matchOp(op, database); });
  }

private:
  llvm::Error loadDatabase(const SliceSpec &spec, MotifDatabase &database) {
    if (!dbPath.empty()) {
      if (!dbSlice.empty()) {
        // db-slice names a fixture slice; silently ignoring it against a JSON
        // database would hide a mistaken invocation.
        return llvm::createStringError(
            llvm::inconvertibleErrorCode(),
            "'db-slice' selects a slice of the compiled-in fixture and cannot "
            "be combined with 'db'");
      }
      return loadJSONDatabase(dbPath, spec, database);
    }
    return loadFixtureDatabase(dbSlice, spec, database);
  }

  /// Annotate `op` with the match result.
  ///
  /// TODO: replace the greedy stand-in with the frontier dynamic program.
  void matchOp(acr::CRNGenericOp op, const MotifDatabase &database) {
    QueryCRN query;
    if (!extractQuery(op, query)) {
      op.emitError("could not read 'stoich_matrix'");
      return signalPassFailure();
    }

    MLIRContext *context = op.getContext();
    llvm::SmallVector<Attribute, 8> speciesNames(
        query.numSpecies, StringAttr::get(context, ""));
    llvm::SmallVector<Attribute, 8> motifIds(query.numReactions,
                                             StringAttr::get(context, ""));
    llvm::SmallVector<unsigned, 4> unmatched;

    for (int64_t r = 0; r < query.numReactions; ++r) {
      unsigned reaction = static_cast<unsigned>(r);
      std::optional<MotifIdx> matched = matchReaction(query, reaction, database);
      if (!matched) {
        unmatched.push_back(reaction);
        continue;
      }

      const Motif &motif = database.getMotif(*matched);
      motifIds[r] = StringAttr::get(context, motif.id);

      // Bind this reaction's species positionally against the motif. The real
      // algorithm searches over injective extensions; the stand-in takes the
      // motif's own order, which is enough to exercise the annotation path.
      llvm::SmallVector<unsigned, 8> supp = query.support(reaction);
      for (auto [i, column] : llvm::enumerate(supp)) {
        if (i >= motif.species.size())
          break;
        auto existing = llvm::cast<StringAttr>(speciesNames[column]);
        if (!existing.getValue().empty())
          continue;
        speciesNames[column] = StringAttr::get(
            context, database.getSpeciesName(motif.species[i]));
      }
    }

    op->setAttr("species_names", ArrayAttr::get(context, speciesNames));
    op->setAttr("motif_ids", ArrayAttr::get(context, motifIds));

    if (!unmatched.empty()) {
      // Partial coverage is a warning, not an error: the paper's flow expects
      // downstream passes to handle partially instantiated CRNs.
      std::string list;
      llvm::raw_string_ostream os(list);
      llvm::interleaveComma(unmatched, os);
      op.emitWarning() << "unmatched reactions: " << os.str();
    }
  }

  /// Pick a motif for one reaction, ignoring cross-reaction consistency.
  std::optional<MotifIdx> matchReaction(const QueryCRN &query, unsigned reaction,
                                        const MotifDatabase &database) {
    llvm::SmallVector<int32_t, 8> coefficients;
    for (int32_t v : query.rows[reaction])
      if (v != 0)
        coefficients.push_back(v);
    if (coefficients.empty())
      return std::nullopt;

    // makeSignature() interns into the database's arena, so build the key
    // through a local database-independent sort instead.
    llvm::sort(coefficients);
    SignatureKey key(coefficients);

    llvm::ArrayRef<MotifIdx> candidates =
        strictSignature ? database.lookupExact(key)
                        : database.lookupSupersets(key);
    if (candidates.empty())
      return std::nullopt;

    // Candidates are in ascending motif-id order, so taking the first is
    // deterministic.
    return candidates.front();
  }
};

} // namespace
} // namespace mlir::crn2db
