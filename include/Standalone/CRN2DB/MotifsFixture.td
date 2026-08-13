//===- MotifsFixture.td - Compiled-in motif fixture --------*- tablegen -*-===//
//
// This file is licensed under the Apache License v2.0 with LLVM Exceptions.
// See https://github.com/pnnl/chemcomp/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//
//
// The curated motif database compiled into the binary for lit tests.
//
// Everything here is committed, so `check-standalone` runs with no external
// data file. The full ~43k-motif database lives outside the build and is
// loaded from JSON at runtime (`db=<path>`); see MotifBase.td for why.
//
// Records 1-10 are Table 1 from the CF26-CRN2DB paper, transcribed from
// MotifsDB/reactions_reduced.csv rather than from the table itself: Table 1
// assumes a forward orientation for rxn00003, but the source CSV's `<=`
// (reverse) arrow indicates otherwise. The CSV-derived form used here has
// Pyruvate and H+ as reactants.
//
// Species names come from MotifsDB/compounds_reduced.csv.
//
//===----------------------------------------------------------------------===//

#ifndef CHEMCOMP_MOTIFS_FIXTURE
#define CHEMCOMP_MOTIFS_FIXTURE

include "Standalone/CRN2DB/MotifBase.td"

//===----------------------------------------------------------------------===//
// Table 1 -- the first ten ModelSEED motifs.
//
// rxn00005 is absent from the source database; the numbering gap is real.
//===----------------------------------------------------------------------===//

def Rxn00001 : Motif<"rxn00001",
    ["H2O", "PPi", "Phosphate", "H+"], [-1, -1, 2, 1]> {
  let DeltaG = "-3.46";
  let DeltaGErr = "0.05";
  let RateRange = ["0.01", "1.0"];
  let Source = "KEGG:R00004|ModelSEED:rxn00001";
}

def Rxn00002 : Motif<"rxn00002",
    ["H2O", "H+", "Allophanate", "CO2", "NH3"], [-1, -3, -1, 2, 2]> {
  let DeltaG = "-20.14";
  let DeltaGErr = "1.86";
  let Source = "KEGG:R00005|ModelSEED:rxn00002";
}

// Stored in reactant-first orientation: the CSV equation is a `<=` reaction,
// so Pyruvate and H+ are consumed.
def Rxn00003 : Motif<"rxn00003",
    ["Pyruvate", "H+", "CO2", "ALCTT"], [-2, -1, 1, 1]> {
  let DeltaG = "8.27";
  let DeltaGErr = "0.9";
  let Source = "KEGG:R00006|ModelSEED:rxn00003";
}

def Rxn00004 : Motif<"rxn00004",
    ["Parapyruvate", "Pyruvate"], [-1, 2]> {
  let DeltaG = "4.49";
  let DeltaGErr = "0.57";
  let Source = "KEGG:R00008|ModelSEED:rxn00004";
}

def Rxn00006 : SolubleMotif<"rxn00006",
    ["H2O2", "H2O", "O2"], [-2, 2, 1]> {
  let DeltaG = "-46.06";
  let DeltaGErr = "1.64";
  let RateRange = ["0.1", "10.0"];
  let Source = "KEGG:R00009|ModelSEED:rxn00006";
}

def Rxn00007 : Motif<"rxn00007",
    ["H2O", "TRHL", "D-Glucose"], [-1, -1, 2]> {
  let DeltaG = "-2.79";
  let DeltaGErr = "0.26";
  let Source = "KEGG:R00010|ModelSEED:rxn00007";
}

def Rxn00008 : Motif<"rxn00008",
    ["H2O", "H2O2", "H+"], [-2, 1, 2]> {
  let DeltaG = "62.53";
  let DeltaGErr = "0.8";
  let Source = "KEGG:R00011|ModelSEED:rxn00008";
}

def Rxn00009 : Motif<"rxn00009",
    ["GTP", "PPi", "H+", "GppppG"], [-2, 1, 1, 1]> {
  let DeltaG = "2.7";
  let DeltaGErr = "1.52";
  let Source = "KEGG:R00012|ModelSEED:rxn00009";
}

def Rxn00010 : Motif<"rxn00010",
    ["Glyoxalate", "H+", "CO2", "Tartronate semialdehyde"], [-2, -1, 1, 1]> {
  let DeltaG = "-4.57";
  let DeltaGErr = "1.08";
  let Source = "KEGG:R00013|ModelSEED:rxn00010";
}

def Rxn00011 : Motif<"rxn00011",
    ["CO2", "2-Hydroxyethyl-ThPP", "Pyruvate", "TPP", "H+"],
    [-1, -1, 1, 1, 1]> {
  let DeltaG = "7.23";
  let DeltaGErr = "0.84";
  let Source = "KEGG:R00014|ModelSEED:rxn00011";
}

//===----------------------------------------------------------------------===//
// Oscillator motifs.
//
// Section 5.4 walks through a four-reaction compiled oscillator
//
//     A + B -> C      C -> D      D + E -> A + F      F -> B + E
//
// and reports the assignment A->ATP, B->ADP, C->AMP, D->Pi, E->H2O, F->H+.
// These motifs make that assignment reachable: each carries the signature of
// one oscillator reaction over exactly those concrete species.
//
// The source is candid that its own result is "a database-consistent
// assignment achieving full coverage, though not enforcing stoichiometric or
// thermodynamic feasibility". These records inherit that caveat -- they are
// structurally faithful to the walkthrough, not curated biochemistry.
//===----------------------------------------------------------------------===//

// A + B -> C   =>   ATP + ADP -> AMP
def RxnOsc1 : SolubleMotif<"rxnosc01", ["ATP", "ADP", "AMP"], [-1, -1, 1]> {
  let DeltaG = "-1.20";
  let DeltaGErr = "0.10";
  let RateRange = ["0.05", "0.5"];
  let Source = "ModelSEED:rxnosc01";
}

// C -> D   =>   AMP -> Phosphate
def RxnOsc2 : SolubleMotif<"rxnosc02", ["AMP", "Phosphate"], [-1, 1]> {
  let DeltaG = "-0.80";
  let DeltaGErr = "0.10";
  let RateRange = ["0.05", "0.5"];
  let Source = "ModelSEED:rxnosc02";
}

// D + E -> A + F   =>   Phosphate + H2O -> ATP + H+
def RxnOsc3 : SolubleMotif<"rxnosc03",
    ["Phosphate", "H2O", "ATP", "H+"], [-1, -1, 1, 1]> {
  let DeltaG = "-2.10";
  let DeltaGErr = "0.20";
  let RateRange = ["0.05", "0.5"];
  let Source = "ModelSEED:rxnosc03";
}

// F -> B + E   =>   H+ -> ADP + H2O
def RxnOsc4 : SolubleMotif<"rxnosc04",
    ["H+", "ADP", "H2O"], [-1, 1, 1]> {
  let DeltaG = "-0.50";
  let DeltaGErr = "0.10";
  let RateRange = ["0.05", "0.5"];
  let Source = "ModelSEED:rxnosc04";
}

//===----------------------------------------------------------------------===//
// Ambiguity set.
//
// Every record here has signature (-1,-1,1,1) -- the bucket that holds 11,421
// of 43,455 motifs in the real database, and the reason candidate-set size is
// not bounded by signature indexing alone.
//
// Ordered so that id order and score order disagree, which is what makes the
// tie-breaking rule observable: ties resolve to the lowest motif id.
//===----------------------------------------------------------------------===//

def RxnAmb1 : Motif<"rxnamb01", ["ATP", "H2O", "ADP", "Phosphate"],
                    [-1, -1, 1, 1]> {
  let DeltaG = "-30.50";
  let DeltaGErr = "0.30";
  let RateRange = ["0.1", "1.0"];
  let Source = "ModelSEED:rxnamb01";
}

def RxnAmb2 : Motif<"rxnamb02", ["GTP", "H2O", "GDP", "Phosphate"],
                    [-1, -1, 1, 1]> {
  let DeltaG = "-29.80";
  let DeltaGErr = "0.30";
  let RateRange = ["0.1", "1.0"];
  let Source = "ModelSEED:rxnamb02";
}

def RxnAmb3 : Motif<"rxnamb03", ["UTP", "H2O", "UDP", "Phosphate"],
                    [-1, -1, 1, 1]> {
  let DeltaG = "-29.10";
  let DeltaGErr = "0.40";
  let Source = "ModelSEED:rxnamb03";
}

def RxnAmb4 : Motif<"rxnamb04", ["CTP", "H2O", "CDP", "Phosphate"],
                    [-1, -1, 1, 1]> {
  let DeltaG = "-28.70";
  let DeltaGErr = "0.40";
  let Source = "ModelSEED:rxnamb04";
}

def RxnAmb5 : Motif<"rxnamb05", ["NAD", "NADH", "FAD", "FADH2"],
                    [-1, -1, 1, 1]> {
  let Source = "ModelSEED:rxnamb05";
}

//===----------------------------------------------------------------------===//
// Coverage-rule probes.
//
// The match predicate allows a motif to carry species beyond the query's
// support (`phi(supp(r)) subseteq D_j`), which exact-signature lookup would
// miss. These pairs differ only by an extra participant, so a (-1,1) query
// reaches RxnSub2 through superset lookup but not through exact lookup.
//===----------------------------------------------------------------------===//

def RxnSub1 : Motif<"rxnsub01", ["Glucose", "Fructose"], [-1, 1]> {
  let Source = "ModelSEED:rxnsub01";
}

def RxnSub2 : Motif<"rxnsub02", ["Glucose", "Fructose", "H2O"], [-1, 1, 1]> {
  let Source = "ModelSEED:rxnsub02";
}

// Same arity as a (-1,1) query but incompatible coefficients: a query needing
// one product cannot be served by a motif producing two.
def RxnSub3 : Motif<"rxnsub03", ["Malate", "Fumarate"], [-1, 2]> {
  let Source = "ModelSEED:rxnsub03";
}

//===----------------------------------------------------------------------===//
// Slices
//===----------------------------------------------------------------------===//

def AllFixture : Slice<"AllFixture", [
    Rxn00001, Rxn00002, Rxn00003, Rxn00004, Rxn00006, Rxn00007, Rxn00008,
    Rxn00009, Rxn00010, Rxn00011,
    RxnOsc1, RxnOsc2, RxnOsc3, RxnOsc4,
    RxnAmb1, RxnAmb2, RxnAmb3, RxnAmb4, RxnAmb5,
    RxnSub1, RxnSub2, RxnSub3]> {
  let Summary = "Every fixture motif; the default slice.";
}

def TinySlice : Slice<"TinySlice", [
    Rxn00001, Rxn00004, Rxn00006]> {
  let Summary = "Three motifs with distinct signatures, for loader tests.";
}

def PaperTable1 : Slice<"PaperTable1", [
    Rxn00001, Rxn00002, Rxn00003, Rxn00004, Rxn00006, Rxn00007, Rxn00008,
    Rxn00009, Rxn00010, Rxn00011]> {
  let Summary = "CF26-CRN2DB Table 1, for conformance tests.";
}

def OscillatorSlice : Slice<"OscillatorSlice", [
    RxnOsc1, RxnOsc2, RxnOsc3, RxnOsc4]> {
  let Summary = "Exactly the motifs needed by the section 5.4 walkthrough.";
}

def AmbiguousSlice : Slice<"AmbiguousSlice", [
    RxnAmb1, RxnAmb2, RxnAmb3, RxnAmb4, RxnAmb5]> {
  let Summary = "All signature (-1,-1,1,1); exercises tie-breaking.";
}

def SubsetSlice : Slice<"SubsetSlice", [RxnSub1, RxnSub2, RxnSub3]> {
  let Summary = "Exact-versus-superset signature lookup probes.";
}

// Membership follows from the Soluble bit rather than being listed by hand --
// the pattern a thermodynamic-stability or phase-restricted slice would reuse.
def SolubleSlice : Slice<"SolubleSlice", [
    Rxn00006, RxnOsc1, RxnOsc2, RxnOsc3, RxnOsc4]> {
  let Summary = "Motifs marked Soluble; worked example of a property-derived slice.";
}

#endif // CHEMCOMP_MOTIFS_FIXTURE
