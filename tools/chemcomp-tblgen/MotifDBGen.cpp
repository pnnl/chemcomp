//===- MotifDBGen.cpp - Motif database generator ----------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Emits a static C++ table from records defined against MotifBase.td.
//
// The output is a self-contained .inc consumed by MotifDatabaseFixture.cpp; it
// declares no types of its own beyond the local record structs, so it can be
// included exactly once without polluting any header.
//
//===----------------------------------------------------------------------===//

#include "llvm/ADT/STLExtras.h"
#include "llvm/ADT/StringExtras.h"
#include "llvm/ADT/StringRef.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/TableGen/Error.h"
#include "llvm/TableGen/Record.h"
#include "llvm/TableGen/TableGenBackend.h"

#include "MotifDBGen.h"

#include <string>
#include <vector>

using namespace llvm;

namespace {

/// Emit `str` as a C++ string literal, escaping what must be escaped.
/// Species names in the real database contain characters like `"` and `\`
/// rarely, but compound names definitely contain commas and parentheses.
void emitStringLiteral(raw_ostream &os, StringRef str) {
  os << '"';
  for (char c : str) {
    switch (c) {
    case '"':
      os << "\\\"";
      break;
    case '\\':
      os << "\\\\";
      break;
    case '\n':
      os << "\\n";
      break;
    case '\t':
      os << "\\t";
      break;
    default:
      os << c;
    }
  }
  os << '"';
}

/// Parse a decimal string field into a double literal, or emit the
/// quiet-NaN spelling when the field is empty.
///
/// Empty means "unknown", which is a distinct state from any numeric value:
/// 43.1% of ModelSEED rows have no thermodynamic data, and treating that as 0
/// would silently fabricate a measurement.
void emitOptionalDouble(raw_ostream &os, StringRef field, Record *rec,
                        StringRef fieldName) {
  if (field.empty()) {
    os << "kUnknown";
    return;
  }
  double value;
  if (field.getAsDouble(value)) {
    PrintFatalError(rec->getLoc(), "field '" + fieldName.str() +
                                       "' is not a valid decimal number: '" +
                                       field.str() + "'");
  }
  // Emit the source spelling rather than a re-formatted double so the
  // generated file stays diffable and no precision is invented.
  os << field;
}

/// Validate one motif record and emit its table row.
void emitMotif(raw_ostream &os, Record *rec) {
  StringRef id = rec->getValueAsString("MotifId");
  std::vector<StringRef> species = rec->getValueAsListOfStrings("Species");
  std::vector<int64_t> stoich = rec->getValueAsListOfInts("Stoich");

  if (species.size() != stoich.size()) {
    PrintFatalError(rec->getLoc(),
                    "motif '" + id.str() + "' has " + Twine(species.size()) +
                        " species but " + Twine(stoich.size()) +
                        " stoichiometric coefficients; the lists must be "
                        "positionally parallel");
  }
  if (species.empty()) {
    PrintFatalError(rec->getLoc(),
                    "motif '" + id.str() + "' has no species");
  }

  // A species repeated within one motif makes nu_j ambiguous. Catching it here
  // gives a source location; MotifDatabase::addMotif would only return false.
  for (size_t i = 0; i < species.size(); ++i) {
    for (size_t j = i + 1; j < species.size(); ++j) {
      if (species[i] == species[j]) {
        PrintFatalError(rec->getLoc(), "motif '" + id.str() +
                                           "' repeats species '" +
                                           species[i].str() + "'");
      }
    }
  }

  std::vector<StringRef> rateRange = rec->getValueAsListOfStrings("RateRange");
  if (!rateRange.empty() && rateRange.size() != 2) {
    PrintFatalError(rec->getLoc(),
                    "motif '" + id.str() + "' has a RateRange of " +
                        Twine(rateRange.size()) +
                        " entries; expected either none or [lo, hi]");
  }

  os << "  {";
  emitStringLiteral(os, id);
  os << ", ";

  // Species names.
  os << "{";
  for (auto [i, s] : enumerate(species)) {
    if (i)
      os << ", ";
    emitStringLiteral(os, s);
  }
  os << "}, ";

  // Coefficients.
  os << "{";
  for (auto [i, v] : enumerate(stoich)) {
    if (i)
      os << ", ";
    os << v;
  }
  os << "}, ";

  emitOptionalDouble(os, rec->getValueAsString("DeltaG"), rec, "DeltaG");
  os << ", ";
  emitOptionalDouble(os, rec->getValueAsString("DeltaGErr"), rec, "DeltaGErr");
  os << ", ";

  // Rate range, as a fixed pair plus a validity flag.
  if (rateRange.empty()) {
    os << "0.0, 0.0, false, ";
  } else {
    emitOptionalDouble(os, rateRange[0], rec, "RateRange[0]");
    os << ", ";
    emitOptionalDouble(os, rateRange[1], rec, "RateRange[1]");
    os << ", true, ";
  }

  emitStringLiteral(os, rec->getValueAsString("Source"));

  // Typed property bits. Adding a column here is the whole cost of extending
  // the schema; nothing downstream of the generated table needs to change.
  unsigned flags = 0;
  if (rec->isSubClassOf("SolubleMotif") && rec->getValueAsBit("Soluble"))
    flags |= 1u << 0;
  os << ", " << flags;

  os << "},\n";
}

} // namespace

bool emitMotifDB(raw_ostream &os, RecordKeeper &records) {
  emitSourceFileHeader("Reaction motif fixture database", os);

  std::vector<Record *> motifs = records.getAllDerivedDefinitions("Motif");
  std::vector<Record *> slices = records.getAllDerivedDefinitions("Slice");

  // Sort by motif id so the generated table is stable regardless of the order
  // definitions appear in the .td, and so it matches the id ordering
  // MotifDatabase::buildIndex() imposes.
  sort(motifs, [](Record *a, Record *b) {
    return a->getValueAsString("MotifId") < b->getValueAsString("MotifId");
  });
  sort(slices, [](Record *a, Record *b) {
    return a->getValueAsString("Name") < b->getValueAsString("Name");
  });

  // Reject duplicate ids: two records claiming one id would make slice
  // membership ambiguous and silently shadow one another.
  for (size_t i = 1; i < motifs.size(); ++i) {
    if (motifs[i]->getValueAsString("MotifId") ==
        motifs[i - 1]->getValueAsString("MotifId")) {
      PrintFatalError(motifs[i]->getLoc(),
                      "duplicate motif id '" +
                          motifs[i]->getValueAsString("MotifId").str() + "'");
    }
  }

  os << "\n#ifdef GET_MOTIF_DB\n";
  os << "#undef GET_MOTIF_DB\n\n";

  os << "// Number of motifs and slices in this fixture.\n";
  os << "static constexpr unsigned kNumFixtureMotifs = " << motifs.size()
     << ";\n";
  os << "static constexpr unsigned kNumFixtureSlices = " << slices.size()
     << ";\n\n";

  os << "// Property bits, matching MotifBase.td's typed subclasses.\n";
  os << "static constexpr unsigned kMotifFlagSoluble = 1u << 0;\n\n";

  os << "namespace {\n";
  os << "constexpr double kUnknown = "
        "std::numeric_limits<double>::quiet_NaN();\n\n";

  os << "struct FixtureMotif {\n";
  os << "  const char *id;\n";
  os << "  std::initializer_list<const char *> species;\n";
  os << "  std::initializer_list<int32_t> stoich;\n";
  os << "  double deltaG;\n";
  os << "  double deltaGErr;\n";
  os << "  double rateLo;\n";
  os << "  double rateHi;\n";
  os << "  bool hasRate;\n";
  os << "  const char *source;\n";
  os << "  unsigned flags;\n";
  os << "};\n\n";

  os << "struct FixtureSlice {\n";
  os << "  const char *name;\n";
  os << "  std::initializer_list<const char *> motifIds;\n";
  os << "};\n\n";

  os << "// Sorted by motif id.\n";
  os << "const FixtureMotif kFixtureMotifs[] = {\n";
  for (Record *rec : motifs)
    emitMotif(os, rec);
  os << "};\n\n";

  os << "// Sorted by slice name; each membership list is sorted by motif id.\n";
  os << "const FixtureSlice kFixtureSlices[] = {\n";
  for (Record *slice : slices) {
    std::vector<Record *> members = slice->getValueAsListOfDefs("Motifs");
    std::vector<StringRef> ids;
    ids.reserve(members.size());
    for (Record *m : members)
      ids.push_back(m->getValueAsString("MotifId"));
    sort(ids);
    ids.erase(std::unique(ids.begin(), ids.end()), ids.end());

    StringRef summary = slice->getValueAsString("Summary");
    if (!summary.empty())
      os << "  // " << summary << "\n";
    os << "  {";
    emitStringLiteral(os, slice->getValueAsString("Name"));
    os << ", {";
    for (auto [i, id] : enumerate(ids)) {
      if (i)
        os << ", ";
      emitStringLiteral(os, id);
    }
    os << "}},\n";
  }
  os << "};\n";
  os << "} // namespace\n\n";
  os << "#endif // GET_MOTIF_DB\n";

  return false;
}
