// Define reaction 1: 2H2 + O2 -> 2H2O
// Reaction 1 adjusts concentrations of H2, O2, and H2O
acr.reaction @reaction_1(%h2: !acr.species, 
                         %o2: !acr.species, 
                         %h2o :!acr.species, 
                         %rate_constant: f32) 
    -> (!acr.species, !acr.species, !acr.species) {
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  %h2_stoich = acr.stoich_coeff -2 : f32
  %o2_stoich = acr.stoich_coeff -1 : f32
  %h2o_stoich = acr.stoich_coeff 2 : f32

  %h2_o, %o2_o, %h2o_o = acr.react(%h2, %o2, %h2o, 
    %h2_stoich, %o2_stoich, %h2o_stoich, 
    %rate_constant) : (!acr.species, !acr.species, !acr.species, 
                       f32, f32, f32, 
                       f32) -> (!acr.species, !acr.species, !acr.species)
  return %h2_o, %o2_o, %h2o_o
}


// Declare a reaction that was defined elsewhere
acr.reaction private @reaction_1(%h2: !acr.species, 
                                 %o2: !acr.species, 
                                 %h2o :!acr.species, 
                                 %rate_constant: f32) 
    -> (!acr.species, !acr.species, !acr.species)

acr.reaction @reaction_2(%co: !acr.species, %o2: !acr.species, 
                         %co2: !acr.species, %rate_2: f32) 
    -> (!acr.species, !acr.species, !acr.species) {
  // Define reaction 2: CO + O2 -> CO2
  // This adjusts concentrations of CO, O2, and CO2
  %co_stoich = acr.stoich_coeff -1 : f32
  %o2_stoich = acr.stoich_coeff -1 : f32
  %co2_stoich = acr.stoich_coeff 1 : f32

  // React with the provided rate constant, adjusting CO, O2, and CO2
  %co_o, %o2_o, %co2_o = acr.react(%co, %o2, %co2, 
                                   %co_stoich, %o2_stoich, %co2_stoich, 
                                   %rate_2) 
    : (!acr.species, !acr.species, !acr.species, f32, f32, f32, f32) 
    -> (!acr.species, !acr.species, !acr.species)
  return %co_o, %o2_o, %co2_o
}

acr.crn @reaction_network(%h2: !acr.species, %o2: !acr.species, 
                          %co: !acr.species, %h2o: !acr.species, 
                          %co2: !acr.species, 
                          %rate_1: f32, %rate_2: f32) 
    -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species) {
  // Instantiate reaction 1: 2H2 + O2 -> 2H2O
  // Reaction 1 adjusts concentrations of H2, O2, and H2O
  %h2_o, %o2_o, %h2o_o = acr.reaction_instance @reaction_1(%h2, %o2, 
                                                           %h2o, %rate_1) 
    : (!acr.species, !acr.species, !acr.species, f32) 
    -> (!acr.species, !acr.species, !acr.species)

  // Instantiate reaction 2 within the network
  %co_o, %o2_o_final, %co2_o = acr.reaction_instance @reaction_2(%co, %o2_o, 
                                                                 %co2, 
                                                                 %rate_2) 
    : (!acr.species, !acr.species, !acr.species, f32) 
    -> (!acr.species, !acr.species, !acr.species)

  // Return updated concentrations for all molecules: H2, O2, H2O, CO, and CO2
  return %h2_o, %o2_o_final, %h2o_o, %co_o, %co2_o
}


// Declare, Declare, and use

// Declare reactions that were defined elsewhere
acr.reaction private @reaction_1( %h2: !acr.species, %o2: !acr.species, 
    %h2o :!acr.species, %rate_constant: f32) 
    -> (!acr.species, !acr.species, !acr.species)

acr.reaction @reaction_2(%co: !acr.species, %o2: !acr.species, 
                         %co2: !acr.species, %rate_2: f32) 
    -> (!acr.species, !acr.species, !acr.species) 


// Use reactions in a network
acr.crn @reaction_network(%h2: !acr.species, %o2: !acr.species, 
                          %co: !acr.species, %h2o: !acr.species, 
                          %co2: !acr.species, 
                          %rate_1: f32, %rate_2: f32) 
    -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species) {
  // Instantiate reaction 1: 2H2 + O2 -> 2H2O
  // Reaction 1 adjusts concentrations of H2, O2, and H2O
  %h2_o, %o2_o, %h2o_o = acr.reaction_instance @reaction_1(%h2, %o2, 
                                                           %h2o, %rate_1) 
    : (!acr.species, !acr.species, !acr.species, f32) 
    -> (!acr.species, !acr.species, !acr.species)

  // Instantiate reaction 2 within the network
  %co_o, %o2_o_final, %co2_o = acr.reaction_instance @reaction_2(%co, %o2_o, 
                                                                 %co2, 
                                                                 %rate_2) 
    : (!acr.species, !acr.species, !acr.species, f32) 
    -> (!acr.species, !acr.species, !acr.species)

  // Return updated concentrations for all molecules: H2, O2, H2O, CO, and CO2
  return %h2_o, %o2_o_final, %h2o_o, %co_o, %co2_o
}


// ----------------
// Define reactions

acr.reaction @reaction_1(%species_vec: vector<3x!acr.species>, 
                         %rate_constant: f32) 
    -> (vector<3x!acr.species>) {
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  %stoich_coeffs = vector.from_elements [-2, -1, 2] : vector<3xf32>

  // Perform the reaction using the species vector and stoichiometry coefficients
  %species_vec_out = acr.react(%species_vec, %stoich_coeffs, %rate_constant) 
    : (vector<3x!acr.species>, vector<3xf32>, f32) 
    -> vector<3x!acr.species>

  return %species_vec_out
}

acr.reaction @reaction_2(%subvector: vector<3x!acr.species>, 
                         %rate_constant: f32) 
    -> (vector<3x!acr.species>) {
  // Define the stoichiometry: CO + O2 -> CO2
  %stoich_coeffs = vector.from_elements [-1, -1, 1] : vector<3xf32>

  // React with the provided rate constant
  %subvector_out = acr.react(%subvector, %stoich_coeffs, %rate_constant) 
    : (vector<3x!acr.species>, vector<3xf32>, f32) 
    -> vector<3x!acr.species>

  return %subvector_out
}

// ------------------------
// Define the CRN operation
acr.crn @reaction_network(%species_vec: vector<5x!acr.species>, 
                          %rate_vec: vector<2xf32>) 
    -> (vector<5x!acr.species>) {
  // For these 5 species: H2, H2O, O2, CO, CO2 
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  //    and the stoichiometry: CO + O2 -> CO2

  // Select the subvector for reaction 1
  %subvector_r1 = vector.gather %species_vec, [0, 1, 2] 
    : vector<5x!acr.species>, vector<3xindex>
  // Instantiate reaction 1 with the rate constant as an input
  %subvector_r1_out = acr.reaction_instance @reaction_1(%subvector_r1, 
                                                        %rate_vec[0]) 
    : (vector<3x!acr.species>, f32) 
    -> vector<3x!acr.species>
  // Update the original species vector with the output subvector
  %species_vec_out1 = vector.scatter %species_vec, %subvector_r1_out, [0, 1, 2] 
    : vector<5x!acr.species>, vector<3x!acr.species>, vector<3xindex>
  
  // Select the subvector for reaction 2
  %subvector_r2 = vector.gather %species_vec_out1, [3, 1, 4] 
    : vector<5x!acr.species>, vector<3xindex>
  // Instantiate reaction 2 with the rate constant as an input
  %subvector_r2_out = acr.reaction_instance @reaction_2(%subvector_r2, 
                                                        %rate_vec[1]) 
    : (vector<3x!acr.species>, f32) 
    -> vector<3x!acr.species>
  // Update the original species vector with the output subvector
  %species_vec_out = vector.scatter %species_vec_out1, %subvector_r2_out, 
                                    [3, 1, 4] 
    : vector<5x!acr.species>, vector<3x!acr.species>, vector<3xindex>

  return %species_vec_out
}

// A generic operation that can represent any CRN
acr.crn_generic {stoich_matrix = 
[
    // With the following order of species:
    // H2, H2O, O2, CO, CO2 
    [-2,-1, 2, 0, 0], //2H2 + O2 -> 2H2O
    [ 0,-1, 0,-1, 1]  //CO + O2 -> CO2
] : tensor<5x2x>
 } (%species_vec: vector<5x!acr.species>, %rate_vec: vector<2xf32>) 
    -> (vector<5x!acr.species>)