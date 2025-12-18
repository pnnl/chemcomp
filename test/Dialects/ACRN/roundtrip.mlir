// // Define the custom CRN dialect
// dialect CRNDialect {
//   // Define a custom operation for a chemical reaction
//   operation chemreaction(%reactants: vector<?xindex>, %products: vector<?xindex>, %rate: f64) {
//     // Implementation details for a single reaction can be added here
//     // For simplicity, let's assume a simple representation
//     // May want to include reaction rates, constants, etc.
//   }

//   // Define a custom operation for a chemical network
//   // Syntax: chemicalnetwork %reactions, %initial_state, %time, %result
//   // Example: chemicalnetwork [%reaction1, %reaction2], [1, 0], 10.0, %final_state
//   // Represents a network with two reactions and initial state [1, 0]
//   // evolving over 10 units of time
//   operation chemicalnetwork(%reactions: vector<?>, %initial_state: vector<?xindex>,
//                             %time: f64, %result: vector<?xindex>) {
//     // Implementation details for network evolution can be added here
//     // Connect reactions, handle time evolution, etc.
//   }

//   // Define a custom type for representing the concentration of chemical species
//   // For simplicity, using a scalar type, but might want to define a more complex type
//   type concentration = i32
// }


// Example usage of the custom dialect
func.func @simulate_reaction() {
  // Define the reactants and products concentrations
//   2A + B -> 3C
  %reactants = arith.constant dense<[2, 1]> : vector<2xi32>
  %products = arith.constant dense<[3]> : vector<1xi32>
  %rate = arith.constant 1.0 : f64

  // Use the custom chemreaction operation
  %result = acr.chemreaction %reactants, %products %rate: (vector<2xi32>, vector<1xi32>, f64) -> ()

  // ...
  return
}

// Example usage of the custom CRN dialect
func.func @simulate_chemical_network() {
  // Define the reactions
  %rate = arith.constant 1.0 : f64
  %reaction1 = acr.chemreaction [2, 1], [3], %rate : vector<2xi32>, vector<1xi32>, f64
  %reaction2 = acr.chemreaction [1], [2, 1], %rate : vector<1xi32>, vector<2xi32>, f64

  // Define the initial state
  %initial_state = arith.constant dense<[1, 0]> : vector<2xi32>

  // Define the time for network evolution
  %time = arith.constant 10.0 : f64

  // Use the custom chemicalnetwork operation
  %final_state = acrn.network [%reaction1, %reaction2], %initial_state, %time : (vector<?>, vector<?xindex>, f64) -> vector<?xindex>

  // The result %final_state represents the state of the network after evolution
  // The rest of the simulation logic can be added here
  // ...

  return
}


// Example usage of the custom CRN dialect
func.func @simulate_chemical_network() {
  // Define the reactions
  %rate = arith.constant 1.0 : f64
  %r1 = acr.chemreaction [1,0,0], [0,1,0], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>
  %r2 = acr.chemreaction [0,1,0], [0,0,1], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>
  %r3 = acr.chemreaction [0,0,1], [1,0,0], %rate : 
    (vector<3xi32>, vector<3xi32>, f64) -> cr<3>

  // Define the initial state
  %initial_state = arith.constant dense<[1.0, 0.13, 1.86]> : vector<3xf32>

  // Define Abstract CRN
  %crn1 = acrn.network [%r1, %r2, %r3], %initial_state : 
    (cr<3>, cr<3>, cr<3>, vector<3xf32>, f64) -> crn<3x3>

  // The result %final_state represents the state of the network after evolution
  // The rest of the simulation logic can be added here
  // ...

  return
}