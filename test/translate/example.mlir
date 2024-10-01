// RUN: standalone-translate --split-input-file -mlir-to-json %s | FileCheck %s 

//             // ODEs:
//             // 400dy/dt - cos(t/400) = 0
//
//             // Chemical Reactions:
//             // A + B -> C
//             // C -> D
//             {
//                 // Reactions are written in a list.
// CHECK-LABEL:    "reactions": [
// CHECK-NEXT:         {
//                         // // A + B -> C
// CHECK-NEXT:             "species"   : ["A", "B", "C"],
// CHECK-NEXT:             "reactants" : [1, 1, 0],
// CHECK-NEXT:             "products"  : [0, 0, 1],
// CHECK-NEXT:             "rate": 0.1,
// CHECK-NEXT:             "name": "R1"
// CHECK-NEXT:         },
// CHECK-NEXT:         {
//                         // C -> D
// CHECK-NEXT:             "species"   : ["C", "D"],
// CHECK-NEXT:             "reactants" : [1, 0],
// CHECK-NEXT:             "products"  : [0, 1],
// CHECK-NEXT:             "rate": 0.05,
// CHECK-NEXT:             "name": "R2"
// CHECK-NEXT:         }
// CHECK-NEXT:     ],
//             
//                 // Initial populations
// CHECK-NEXT:     "populations": {
// CHECK-NEXT:         "A": 1,
// CHECK-NEXT:         "B": 1,
// CHECK-NEXT:         "C": 0,
// CHECK-NEXT:         "D": 0
// CHECK-NEXT:     },
//             
//                 // List of species to track on the simulator, 
//                 // this is the concentration that should implement our ODE
// CHECK-NEXT:     "tracked_species": ["D"]
// CHECK-NEXT: }


// -----

//             // Chemical Reactions:
//             // A + B -> 2B
//             // B <-> A
//             {
//                 // Reactions are written in a list.
// CHECK-LABEL:    "reactions": [
// CHECK-NEXT:         {
//                         // A + B -> 2B
// CHECK-NEXT:             "species"   : ["A", "B"],
// CHECK-NEXT:             "reactants" : [1, 1],
// CHECK-NEXT:             "products"  : [0, 2],
// CHECK-NEXT:             "rate": 2,
// CHECK-NEXT:             "name": "R1"
// CHECK-NEXT:         },
// CHECK-NEXT:         {
//                         // B -> A (fwd)
// CHECK-NEXT:             "species"   : ["A", "B"],
// CHECK-NEXT:             "reactants" : [0, 1],
// CHECK-NEXT:             "products"  : [1, 0],
// CHECK-NEXT:             "rate": 2,
// CHECK-NEXT:             "name": "R2"
// CHECK-NEXT:         },
// CHECK-NEXT:         {
//                         // A -> B (bwd)
// CHECK-NEXT:             "species"   : ["A", "B"],
// CHECK-NEXT:             "reactants" : [1, 0],
// CHECK-NEXT:             "products"  : [0, 1],
// CHECK-NEXT:             "rate": 0.1,
// CHECK-NEXT:             "name": "R3"
// CHECK-NEXT:         }
// CHECK-NEXT:     ],
//             
//                 // Initial populations
// CHECK-NEXT:     "populations": {
// CHECK-NEXT:         "A": 10,
// CHECK-NEXT:         "B": 5
// CHECK-NEXT:     },
//             
//                 // List of species to track on the simulator, 
//                 // this is the concentration that should implement our ODE
// CHECK-NEXT:     "tracked_species": ["B"]
// CHECK-NEXT: }