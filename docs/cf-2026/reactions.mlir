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
  acr.return %h2_o, %o2_o, %h2o_o : !acr.species, !acr.species, !acr.species
}


// Declare a reaction that was defined elsewhere
acr.reaction private @reaction_1(!acr.species, !acr.species, !acr.species, f32)
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
  acr.return %co_o, %o2_o, %co2_o : !acr.species, !acr.species, !acr.species
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
  acr.return %h2_o, %o2_o_final, %h2o_o, %co_o, %co2_o
      : !acr.species, !acr.species, !acr.species, !acr.species, !acr.species
}


// Declare, Declare, and use

// Declare reactions that were defined elsewhere
acr.reaction private @reaction_1(!acr.species, !acr.species, !acr.species, f32)
    -> (!acr.species, !acr.species, !acr.species)

acr.reaction private @reaction_2(!acr.species, !acr.species, !acr.species, f32)
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
  acr.return %h2_o, %o2_o_final, %h2o_o, %co_o, %co2_o
      : !acr.species, !acr.species, !acr.species, !acr.species, !acr.species
}


// ----------------
// Tensor form: Define reactions using tensor<Nx!acr.species>

acr.reaction @reaction_1_tensor(%species_vec: tensor<3x!acr.species>,
                                %rate_constant: f32)
    -> tensor<3x!acr.species> {
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  %stoich_coeffs = arith.constant dense<[-2.0, -1.0, 2.0]> : tensor<3xf32>

  // Perform the reaction using the species tensor and stoichiometry coefficients
  %species_vec_out = acr.react(%species_vec, %stoich_coeffs, %rate_constant)
    : (tensor<3x!acr.species>, tensor<3xf32>, f32)
    -> tensor<3x!acr.species>

  acr.return %species_vec_out : tensor<3x!acr.species>
}

acr.reaction @reaction_2_tensor(%subvector: tensor<3x!acr.species>,
                                %rate_constant: f32)
    -> tensor<3x!acr.species> {
  // Define the stoichiometry: CO + O2 -> CO2
  %stoich_coeffs = arith.constant dense<[-1.0, -1.0, 1.0]> : tensor<3xf32>

  // React with the provided rate constant
  %subvector_out = acr.react(%subvector, %stoich_coeffs, %rate_constant)
    : (tensor<3x!acr.species>, tensor<3xf32>, f32)
    -> tensor<3x!acr.species>

  acr.return %subvector_out : tensor<3x!acr.species>
}

// ------------------------
// Define the CRN operation using acr.crn_generic with tensor form
acr.crn @reaction_network_tensor(%species_vec: tensor<5x!acr.species>,
                                  %rate_vec: tensor<2xf32>)
    -> tensor<5x!acr.species> {
  // For these 5 species: H2, H2O, O2, CO, CO2
  // Stoichiometric matrix encodes both reactions:
  //   2H2 + O2 -> 2H2O  and  CO + O2 -> CO2
  %species_vec_out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2, 0, 0], [0, -1, 0, -1, 1]]> : tensor<2x5xi32>}
      (%species_vec, %rate_vec)
          : (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>

  acr.return %species_vec_out : tensor<5x!acr.species>
}

// A generic operation that can represent any CRN (scalar form)
acr.crn @generic_network_scalar(%s0: !acr.species, %s1: !acr.species,
                                 %s2: !acr.species, %s3: !acr.species,
                                 %s4: !acr.species,
                                 %r0: f32, %r1: f32)
    -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species) {
  // With the following order of species: H2, H2O, O2, CO, CO2
  %out0, %out1, %out2, %out3, %out4 = acr.crn_generic
      {stoich_matrix = dense<[
          [-2,-1, 2, 0, 0], // 2H2 + O2 -> 2H2O
          [ 0,-1, 0,-1, 1]  // CO + O2 -> CO2
      ]> : tensor<2x5xi32>}
      (%s0, %s1, %s2, %s3, %s4, %r0, %r1)
          : (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species, f32, f32)
         -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species)
  acr.return %out0, %out1, %out2, %out3, %out4
      : !acr.species, !acr.species, !acr.species, !acr.species, !acr.species
}
