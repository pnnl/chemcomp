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

#include "Standalone/Conversion/CRNToDB/CRNToDB.h"

#include "Standalone/Conversion/CRNToDB/CRNMatcher.h"
#include "Standalone/Conversion/CRNToDB/FrontierDP.h"
#include "Standalone/Conversion/CRNToDB/MotifDatabase.h"
#include "Standalone/Conversion/CRNToDB/MotifDatabaseLoader.h"

#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinAttributes.h"
#include "mlir/IR/BuiltinTypes.h"

#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/raw_ostream.h"

using namespace mlir;
using namespace chemcomp;

namespace mlir {
#define GEN_PASS_DEF_CRNDBMATCH
#include "Standalone/Conversion/Passes.h.inc"

namespace {

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

  query.numReactions = static_cast<unsigned>(type.getDimSize(0));
  query.numSpecies = static_cast<unsigned>(type.getDimSize(1));

  auto values = stoichAttr.getValues<llvm::APInt>();
  auto it = values.begin();
  query.rows.reserve(query.numReactions);
  for (unsigned r = 0; r < query.numReactions; ++r) {
    llvm::SmallVector<int32_t, 8> row;
    row.reserve(query.numSpecies);
    for (unsigned c = 0; c < query.numSpecies; ++c, ++it)
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
    // Note that `max-candidates` is deliberately NOT applied here as a
    // per-signature slice cap. Cutting the database at load time would hide
    // the truncation from preprocessing, which is what emits the remark
    // warning that the result has become heuristic. The cap belongs to
    // candidate selection, where it is observable.
    SliceSpec spec;
    for (const std::string &s : filterSpecies)
      spec.requireAnySpecies.push_back(s);
    spec.sourcePrefix = filterSource;
    spec.maxPerSignature = maxPerSignature;

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

  /// Match `op` against the database and annotate it with the result.
  void matchOp(acr::CRNGenericOp op, const MotifDatabase &database) {
    QueryCRN query;
    if (!extractQuery(op, query)) {
      op.emitError("could not read 'stoich_matrix'");
      return signalPassFailure();
    }

    PreprocessOptions preprocessOptions;
    preprocessOptions.maxCandidates = maxCandidates;
    preprocessOptions.strictSignature = strictSignature;
    Preprocessing prep = preprocess(query, database, preprocessOptions);

    if (dumpPreprocessing) {
      prep.print(llvm::outs(), query, database);
      return;
    }

    // Truncation turns an exact algorithm into a heuristic, so it must never
    // be silent.
    for (auto [reaction, before] : prep.truncated)
      op.emitRemark() << "candidate set for reaction " << reaction
                      << " truncated from " << before << " to "
                      << prep.candidates[reaction].size()
                      << "; the result is heuristic rather than optimal";

    if (prep.arcConsistencyIncomplete)
      op.emitRemark() << "arc-consistency stopped at its iteration budget; "
                         "candidate domains are larger than a full run would "
                         "leave them";

    DPOptions dpOptions;
    dpOptions.weights.coverage = weightCoverage;
    dpOptions.weights.stoich = weightStoich;
    dpOptions.weights.rate = weightRate;
    dpOptions.weights.sigma = rateSigma;
    dpOptions.beamWidth = beamWidth;
    dpOptions.globalInjectivity = globalInjectivity;

    MatchResult match = runFrontierDP(query, database, prep, dpOptions);

    MLIRContext *context = op.getContext();
    llvm::SmallVector<Attribute, 8> speciesNames;
    speciesNames.reserve(query.numSpecies);
    for (const std::optional<SpeciesId> &bound : match.assignment)
      speciesNames.push_back(StringAttr::get(
          context, bound ? database.getSpeciesName(*bound) : llvm::StringRef()));

    llvm::SmallVector<Attribute, 8> motifIds;
    motifIds.reserve(query.numReactions);
    for (const std::optional<MotifIdx> &matched : match.matches)
      motifIds.push_back(StringAttr::get(
          context, matched ? database.getMotif(*matched).id : llvm::StringRef()));

    op->setAttr("species_names", ArrayAttr::get(context, speciesNames));
    op->setAttr("motif_ids", ArrayAttr::get(context, motifIds));

    if (match.beamTruncated)
      op.emitRemark() << "beam search discarded states; the result may be "
                         "suboptimal";

    if (!match.injectivityViolations.empty()) {
      // Frontier-local injectivity is the paper's own state design: a binding
      // is forgotten once its species leaves the frontier, so a retired
      // chemical can be reused. Report rather than hide it, and point at the
      // option that enforces the paper's stated requirement.
      std::string list;
      llvm::raw_string_ostream os(list);
      llvm::interleaveComma(match.injectivityViolations, os,
                            [&](SpeciesId s) { os << database.getSpeciesName(s); });
      op.emitWarning() << "assignment is not globally injective; reused: "
                       << os.str()
                       << " (pass global-injectivity=1 to enforce injectivity)";
    }

    if (!match.unmatched.empty()) {
      // Partial coverage is a warning, not an error: the paper's flow expects
      // downstream passes to handle partially instantiated CRNs.
      std::string list;
      llvm::raw_string_ostream os(list);
      llvm::interleaveComma(match.unmatched, os);
      op.emitWarning() << "unmatched reactions: " << os.str();
    }
  }
};

} // namespace
} // namespace mlir
