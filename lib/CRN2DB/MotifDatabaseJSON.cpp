//===- MotifDatabaseJSON.cpp - JSON motif database loader -------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone/CRN2DB/MotifDatabaseLoader.h"

#include "llvm/ADT/SmallVector.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/JSON.h"
#include "llvm/Support/MemoryBuffer.h"

using namespace chemcomp;

namespace {

llvm::Error makeError(const llvm::Twine &msg) {
  return llvm::createStringError(llvm::inconvertibleErrorCode(), msg);
}

/// Read an optional numeric field.
///
/// Absent, null, and the 1e7 ModelSEED sentinel all mean "no measurement" and
/// map to NaN. A present-but-wrong-type value is an error rather than a
/// silent default, because that indicates a malformed file.
llvm::Error readOptionalDouble(const llvm::json::Object &obj,
                               llvm::StringRef key, unsigned index,
                               double &out) {
  out = std::numeric_limits<double>::quiet_NaN();

  const llvm::json::Value *value = obj.get(key);
  if (!value)
    return llvm::Error::success();
  if (value->kind() == llvm::json::Value::Null)
    return llvm::Error::success();

  std::optional<double> number = value->getAsNumber();
  if (!number)
    return makeError("record " + llvm::Twine(index) + ": field '" + key +
                     "' must be a number or null");

  if (!isUnknownDeltaG(*number))
    out = *number;
  return llvm::Error::success();
}

llvm::Error parseRecords(llvm::StringRef json, llvm::StringRef context,
                         const SliceSpec &spec, MotifDatabase &db) {
  llvm::Expected<llvm::json::Value> parsed = llvm::json::parse(json);
  if (!parsed)
    return makeError(context + ": " + llvm::toString(parsed.takeError()));

  const llvm::json::Array *records = parsed->getAsArray();
  if (!records)
    return makeError(context +
                     ": expected a top-level array of motif objects");

  // Reused across records to keep the 43k-record path from reallocating.
  llvm::SmallVector<llvm::StringRef, 8> speciesNames;
  llvm::SmallVector<int32_t, 8> stoich;
  llvm::SmallVector<double, 2> rateRange;

  for (auto [index, entry] : llvm::enumerate(*records)) {
    const llvm::json::Object *obj = entry.getAsObject();
    if (!obj)
      return makeError(context + ": record " + llvm::Twine(index) +
                       " is not an object");

    std::optional<llvm::StringRef> id = obj->getString("motif_id");
    if (!id)
      return makeError(context + ": record " + llvm::Twine(index) +
                       " has no 'motif_id' string");

    const llvm::json::Array *speciesArray = obj->getArray("species");
    if (!speciesArray)
      return makeError(context + ": motif '" + *id + "' has no 'species' array");

    const llvm::json::Array *stoichArray = obj->getArray("stoich");
    if (!stoichArray)
      return makeError(context + ": motif '" + *id + "' has no 'stoich' array");

    if (speciesArray->size() != stoichArray->size())
      return makeError(context + ": motif '" + *id + "' has " +
                       llvm::Twine(speciesArray->size()) + " species but " +
                       llvm::Twine(stoichArray->size()) +
                       " stoichiometric coefficients; the arrays must be "
                       "positionally parallel");

    speciesNames.clear();
    for (const llvm::json::Value &s : *speciesArray) {
      std::optional<llvm::StringRef> name = s.getAsString();
      if (!name)
        return makeError(context + ": motif '" + *id +
                         "' has a non-string species name");
      speciesNames.push_back(*name);
    }

    stoich.clear();
    for (const llvm::json::Value &c : *stoichArray) {
      std::optional<int64_t> coeff = c.getAsInteger();
      if (!coeff)
        return makeError(context + ": motif '" + *id +
                         "' has a non-integer stoichiometric coefficient; "
                         "fractional coefficients are not representable in "
                         "stoich_matrix");
      stoich.push_back(static_cast<int32_t>(*coeff));
    }

    double deltaG, deltaGErr;
    if (llvm::Error e = readOptionalDouble(*obj, "deltaG", index, deltaG))
      return e;
    if (llvm::Error e = readOptionalDouble(*obj, "deltaGerr", index, deltaGErr))
      return e;

    rateRange.clear();
    if (const llvm::json::Array *rates = obj->getArray("rate_range")) {
      if (rates->size() != 2)
        return makeError(context + ": motif '" + *id + "' has a 'rate_range' of " +
                         llvm::Twine(rates->size()) +
                         " entries; expected [lo, hi]");
      for (const llvm::json::Value &r : *rates) {
        std::optional<double> value = r.getAsNumber();
        if (!value)
          return makeError(context + ": motif '" + *id +
                           "' has a non-numeric 'rate_range' entry");
        rateRange.push_back(*value);
      }
    }

    llvm::StringRef source = obj->getString("source").value_or(llvm::StringRef());

    if (!db.addMotif(*id, speciesNames, stoich, deltaG, deltaGErr, rateRange,
                     source)) {
      // addMotif only rejects a length mismatch (checked above) or a species
      // repeated within the motif, which would make nu_j ambiguous.
      return makeError(context + ": motif '" + *id +
                       "' repeats a species; its stoichiometry is ambiguous");
    }
  }

  db.applySlice(spec);
  db.buildIndex();
  return llvm::Error::success();
}

} // namespace

llvm::Error chemcomp::loadJSONDatabaseFromBuffer(llvm::StringRef json,
                                                 llvm::StringRef context,
                                                 const SliceSpec &spec,
                                                 MotifDatabase &db) {
  return parseRecords(json, context, spec, db);
}

llvm::Error chemcomp::loadJSONDatabase(llvm::StringRef path,
                                       const SliceSpec &spec,
                                       MotifDatabase &db) {
  llvm::ErrorOr<std::unique_ptr<llvm::MemoryBuffer>> buffer =
      llvm::MemoryBuffer::getFile(path);
  if (!buffer)
    return makeError("could not open motif database '" + path +
                     "': " + buffer.getError().message());

  return parseRecords((*buffer)->getBuffer(), path, spec, db);
}
