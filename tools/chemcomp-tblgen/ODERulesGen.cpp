//===- ODERulesGen.cpp - ODE rule table generator ---------------*- C++ -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// Emits the rows of the convert-ode-to-acr rule table from ODERules.td, sorted
// by descending priority. Each row is
//
//   {"name", priority, "A + B -> C", "summary",
//    {"A", "B", "C"}, /*reactants=*/{1, 1, 0}, /*products=*/{0, 0, 1}},
//
//===----------------------------------------------------------------------===//

#include "ODERulesGen.h"

#include "llvm/ADT/StringExtras.h"
#include "llvm/Support/raw_ostream.h"
#include "llvm/TableGen/Error.h"
#include "llvm/TableGen/Record.h"
#include "llvm/TableGen/TableGenBackend.h"

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>

using namespace llvm;

namespace {

struct Rule {
  const Record *record;
  std::string name;
  int64_t priority;
  std::string summary;
  std::vector<std::string> species;
  std::vector<unsigned> reactants;
  std::vector<unsigned> products;
};

std::string complexString(const Rule &rule, const std::vector<unsigned> &c) {
  std::string result;
  for (size_t i = 0; i < rule.species.size(); ++i) {
    if (!c[i])
      continue;
    if (!result.empty())
      result += " + ";
    if (c[i] != 1)
      result += std::to_string(c[i]);
    result += rule.species[i];
  }
  return result.empty() ? "0" : result;
}

Rule readRule(const Record *rec) {
  Rule rule{rec, rec->getValueAsString("Name").str(),
            rec->getValueAsInt("Priority"),
            rec->getValueAsString("Summary").str(), {}, {}, {}};
  std::vector<StringRef> reactants = rec->getValueAsListOfStrings("Reactants");
  std::vector<StringRef> products = rec->getValueAsListOfStrings("Products");
  // Species in order of first appearance, reactants first.
  auto indexOf = [&](StringRef name) {
    auto it = std::find(rule.species.begin(), rule.species.end(), name.str());
    if (it != rule.species.end())
      return static_cast<size_t>(it - rule.species.begin());
    rule.species.push_back(name.str());
    rule.reactants.push_back(0);
    rule.products.push_back(0);
    return rule.species.size() - 1;
  };
  for (StringRef name : reactants)
    ++rule.reactants[indexOf(name)];
  for (StringRef name : products)
    ++rule.products[indexOf(name)];

  if (rule.reactants == rule.products)
    PrintFatalError(rec->getLoc(), "rule '" + rule.name +
                                       "' has the same reactants and products");
  if (rule.species.size() > 4)
    PrintFatalError(rec->getLoc(), "rule '" + rule.name +
                                       "' has more than 4 species; matching "
                                       "enumerates species assignments");
  return rule;
}

void emitList(raw_ostream &os, const std::vector<unsigned> &values) {
  os << "{";
  interleaveComma(values, os);
  os << "}";
}

} // namespace

bool emitODERules(raw_ostream &os, const RecordKeeper &records) {
  emitSourceFileHeader("Composable ODE-to-reaction rules", os, records);

  std::vector<Rule> rules;
  std::set<std::string> names;
  std::map<int64_t, std::string> priorities;
  for (const Record *rec : records.getAllDerivedDefinitions("ODERule")) {
    Rule rule = readRule(rec);
    if (!names.insert(rule.name).second)
      PrintFatalError(rec->getLoc(), "duplicate rule name '" + rule.name + "'");
    auto [it, inserted] = priorities.emplace(rule.priority, rule.name);
    if (!inserted)
      PrintFatalError(rec->getLoc(), "rule '" + rule.name +
                                         "' has the same priority as '" +
                                         it->second +
                                         "'; priorities fix the greedy order "
                                         "and must be distinct");
    rules.push_back(std::move(rule));
  }
  std::sort(rules.begin(), rules.end(), [](const Rule &a, const Rule &b) {
    return a.priority > b.priority;
  });

  for (const Rule &rule : rules) {
    os << "{\"" << rule.name << "\", " << rule.priority << ", \""
       << complexString(rule, rule.reactants) << " -> "
       << complexString(rule, rule.products) << "\", \"" << rule.summary
       << "\",\n {";
    interleaveComma(rule.species, os,
                    [&](const std::string &s) { os << '"' << s << '"'; });
    os << "}, ";
    emitList(os, rule.reactants);
    os << ", ";
    emitList(os, rule.products);
    os << "},\n";
  }
  return false;
}
