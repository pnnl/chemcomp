// The whole point of the pass: matched species names must reach the emulator.
//
// crn-db-match annotates the operation, then the existing JSON export carries
// those names out of MLIR. Before this worked, the export emitted generic
// S0, S1, ... and the concrete chemistry stopped at the IR boundary.

// RUN: chemcomp-opt %s -crn-db-match='db-slice=OscillatorSlice' \
// RUN:   | chemcomp-translate -mlir-to-json | FileCheck %s

// A CRN that was never matched still exports, falling back to generic names.
// RUN: chemcomp-translate -mlir-to-json %s | FileCheck %s --check-prefix=UNMATCHED

acr.crn @oscillator(%sv: tensor<6x!acr.species>, %rv: tensor<4xf32>)
    -> tensor<6x!acr.species> {
  %out = acr.crn_generic
      {stoich_matrix = dense<[[-1, -1,  1,  0,  0,  0],
                              [ 0,  0, -1,  1,  0,  0],
                              [ 1,  0,  0, -1, -1,  1],
                              [ 0,  1,  0,  0,  1, -1]]> : tensor<4x6xi32>}
      (%sv, %rv) : (tensor<6x!acr.species>, tensor<4xf32>)
                -> tensor<6x!acr.species>
  acr.return %out : tensor<6x!acr.species>
}

// CHECK: "species": [
// CHECK-NEXT: "ATP"
// CHECK-NEXT: "ADP"
// CHECK-NEXT: "AMP"
// CHECK-NEXT: "Phosphate"
// CHECK-NEXT: "H2O"
// CHECK-NEXT: "H+"

// Initial populations are keyed by the same names.
// CHECK: "populations": {
// CHECK-NEXT: "ATP"

// UNMATCHED: "species": [
// UNMATCHED-NEXT: "S0"
// UNMATCHED-NEXT: "S1"
