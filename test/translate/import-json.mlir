// RUN: chemcomp-translate -json-to-mlir %S/Inputs/import-test.crn.json | FileCheck %s

// CHECK-LABEL: acr.crn @crn
// CHECK-SAME: (%{{.*}}: tensor<3x!acr.species>, %{{.*}}: tensor<2xf32>)
// CHECK-SAME: -> tensor<3x!acr.species>
// CHECK:        acr.crn_generic
// CHECK-SAME:     stoich_matrix = dense<{{\[}}{{\[}}-2, -1, 2], [0, -1, 1]]>
// CHECK-SAME:     tensor<2x3xi32>
// CHECK-SAME:     tensor<3x!acr.species>, tensor<2xf32>
// CHECK-SAME:     -> tensor<3x!acr.species>
// CHECK:        acr.return
