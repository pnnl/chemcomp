## Introduction to Chemical Reactions and ChemComp

Chemical reactions are processes in which one or more substances, known as reactants, are transformed into different substances, known as products. These transformations are based on the rearrangement of atoms and the making or breaking of chemical bonds, without altering the total number of each type of atom involved. In essence, chemical reactions obey the conservation of mass, charge, and energy. In ChemComp we leverage the chemical reaction properties, described in the following sections, and implement a compilation flow capable of transforming  mathematical equations into chemical reaction networks.

### Overview of Chemical Reactions
At the molecular level, chemical reactions involve interactions between atoms, molecules, or ions that result in the formation or breaking of chemical bonds. The reaction mechanism describes the step-by-step changes, including intermediate compounds and transition states, but the general form of chemical reactions can be written as:

$$ \text{Reactants} \rightarrow \text{Products} $$

Chemical reactions can be classified into several types:
- **Synthesis reactions**: Two or more simple substances combine to form a more complex substance.
    $ A + B \rightarrow AB $
- **Decomposition reactions**: A compound breaks down into simpler compounds or elements.
    $ AB \rightarrow A + B $
- **Single-displacement reactions**: One element replaces another element in a compound.
    $ A + BC \rightarrow AC + B $
- **Double-displacement reactions**: Two compounds exchange ions to form new compounds.
    $ AB + CD \rightarrow AD + CB $
- **Redox reactions**: Transfer of electrons between species, where oxidation and reduction occur.
- **Catalytic reactions**: A catalyst speeds up the reaction by lowering the activation energy without being consumed in the process. It facilitates the conversion of reactants to products.
    $ A + B \xrightarrow{\text{Catalyst}} C + D $

### Stoichiometric Equations
A **stoichiometric equation** represents a chemical reaction, indicating both the reactants and products, along with their relative amounts. For instance:

$ aA + bB \rightarrow cC + dD $

Here, \(A\) and \(B\) are reactants, while \(C\) and \(D\) are products. The coefficients \(a\), \(b\), \(c\), and \(d\) represent the stoichiometric coefficients, which are the ratios in which reactants combine and products form. These coefficients are determined by balancing the chemical equation according to the **law of conservation of mass**—the same number of atoms for each element must exist on both sides of the equation.

#### Key Elements of Stoichiometric Equations:
- **Conservation of atoms**: The total number of each type of atom on the reactant side must equal the total number of each type on the product side.
- **Mole ratios**: Stoichiometric coefficients indicate the molar ratio between reactants and products, which is crucial for predicting the outcome of reactions and determining the amounts of reactants needed or products formed.

Example of a stoichiometric equation:
$$ 2H_2 + O_2 \rightarrow 2H_2O $$
This indicates that two molecules of hydrogen react with one molecule of oxygen to produce two molecules of water.

### Representation of Chemical Reactions in Computation
Representing chemical reactions computationally requires understanding how to encode:
1. **Reactants and products**: Represent molecules with atoms and bonds (e.g., graph representations where vertices are atoms and edges are bonds).
2. **Reaction mechanisms**: These can be represented as transitions between molecular graphs.
3. **Stoichiometry**: Ensure that the coefficients in the equation respect mass and charge balance.
4. **Kinetics and thermodynamics**: Additional information might include rate constants for reaction rates and energy changes.

In computational models, chemical reactions are often represented in the form of **reaction networks** or **graphs**. Each node might represent a chemical species, and each edge represents a reaction. The directionality of the edges is important for showing which species are reactants and which are products.

### Abstract Representation for Reasoning
For the purpose of reasoning about abstract chemical reactions using MLIR (Multi-Level Intermediate Representation), you would likely be dealing with a graph-based or symbolic representation. In this framework, chemical species can be represented as nodes, while the transformations (reactions) are edges connecting the nodes. 

The MLIR dialect you design might abstract:
- **Species identifiers**: Unique symbols or identifiers representing chemical compounds.
- **Reactions**: Directed edges annotated with stoichiometry and possibly additional metadata like rate constants or conditions.
- **Balancing constraints**: Mechanisms to enforce stoichiometric balance based on conservation laws.

### Example Considerations for MLIR Infrastructure
1. **Type system**: You may want to create a type system where basic entities are atoms or functional groups, and molecules are compositions of these entities.
2. **Transformation rules**: Define transformations between molecular graphs according to reaction rules. This could involve pattern-matching mechanisms in the graph.
3. **Constraints**: Enforce mass and charge conservation using equations or symbolic constraints in the IR.
4. **Modularity**: Design the dialect to be extensible for different types of reactions (e.g., equilibrium reactions, redox reactions).

### Design principles for MLIR Infrastructure

When designing MLIR operations for chemical reactions, it is essential to distinguish between elements that change (inputs) and elements that remain constant (constraints). For instance, in a chemical reaction, species concentrations and reaction rates might vary, while the stoichiometry of the reaction remains fixed. This separation allows for more flexible and efficient modeling and computation.

### Example of Abstract Chemical Reaction in ChemComp

In ChemComp, we create a dialect for abstract chemical reactions, called `acr`. Here is an example of how a chemical reaction might be represented in MLIR using the ChemComp dialect, where species concentrations and reaction rates are inputs, and stoichiometry is integrated as a constraint:

```mlir
// Define reaction 1: 2H2 + O2 -> 2H2O
// Reaction 1 adjusts concentrations of H2, O2, and H2O
acr.reaction @reaction_1(%h2: !acr.species, %o2: !acr.species, %h2o :!acr.species, %rate_constant: f32) -> (!acr.species, !acr.species, !acr.species) {
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  %h2_stoich = acr.stoich_coeff -2 : f32
  %o2_stoich = acr.stoich_coeff -1 : f32
  %h2o_stoich = acr.stoich_coeff 2 : f32

  %h2_o, %o2_o, %h2o_o = acr.react(%h2, %o2, %h2o, %h2_stoich, %o2_stoich, %h2o_stoich, %rate_constant) : (!acr.species, !acr.species, !acr.species, f32, f32, f32, f32) -> (!acr.species, !acr.species, !acr.species)
  acr.return %h2_o, %o2_o, %h2o_o : !acr.species, !acr.species, !acr.species
}
```

A more elaborate example is related to a network of multiple reactions, a chemical reaction network.

```mlir

// Declare a reaction that was defined elsewhere
acr.reaction private @reaction_1(!acr.species, !acr.species, !acr.species, f32) -> (!acr.species, !acr.species, !acr.species)

acr.reaction @reaction_2(%co: !acr.species, %o2: !acr.species, %co2: !acr.species, %rate_2: f32) -> (!acr.species, !acr.species, !acr.species) {
  // Define reaction 2: CO + O2 -> CO2
  // This adjusts concentrations of CO, O2, and CO2
  %co_stoich = acr.stoich_coeff -1 : f32
  %o2_stoich = acr.stoich_coeff -1 : f32
  %co2_stoich = acr.stoich_coeff 1 : f32

  // React with the provided rate constant, adjusting CO, O2, and CO2
  %co_o, %o2_o, %co2_o = acr.react(%co, %o2, %co2, %co_stoich, %o2_stoich, %co2_stoich, %rate_2) : (!acr.species, !acr.species, !acr.species, f32, f32, f32, f32) -> (!acr.species, !acr.species, !acr.species)
  acr.return %co_o, %o2_o, %co2_o : !acr.species, !acr.species, !acr.species
}

acr.crn @reaction_network(%h2: !acr.species, %o2: !acr.species, %co: !acr.species, %h2o: !acr.species, %co2: !acr.species, %rate_1: f32, %rate_2: f32) -> (!acr.species, !acr.species, !acr.species, !acr.species, !acr.species) {

  // Instantiate reaction 1: 2H2 + O2 -> 2H2O
  // Reaction 1 adjusts concentrations of H2, O2, and H2O
  %h2_o, %o2_o, %h2o_o = acr.reaction_instance @reaction_1(%h2, %o2, %h2o, %rate_1) : (!acr.species, !acr.species, !acr.species, f32) -> (!acr.species, !acr.species, !acr.species)

  // Instantiate reaction 2 within the network
  %co_o, %o2_o_final, %co2_o = acr.reaction_instance @reaction_2(%co, %o2_o, %co2, %rate_2) : (!acr.species, !acr.species, !acr.species, f32) -> (!acr.species, !acr.species, !acr.species)

  // Return updated concentrations for all molecules: H2, O2, H2O, CO, and CO2
  acr.return %h2_o, %o2_o_final, %h2o_o, %co_o, %co2_o : !acr.species, !acr.species, !acr.species, !acr.species, !acr.species
}
```

In this example, the `acr.crn` operation defines a network of reactions involving hydrogen, oxygen, and carbon monoxide. The network includes two reactions, each with its own rate constant. The network is defined as a composition of these reactions, and the results are returned as the output of the network.


### Alternative representation

Another way to use the `acr.crn` operation is with tensors of molecules and rate constants, which can be useful for simulating large reaction networks:

```mlir

// ----------------
// Define reactions

acr.reaction @reaction_1(%species_vec: tensor<3x!acr.species>, %rate_constant: f32) -> tensor<3x!acr.species> {
  // Define the stoichiometry: 2H2 + O2 -> 2H2O
  %stoich_coeffs = arith.constant dense<[-2.0, -1.0, 2.0]> : tensor<3xf32>

  // Perform the reaction using the species tensor and stoichiometry coefficients
  %species_vec_out = acr.react(%species_vec, %stoich_coeffs, %rate_constant) : (tensor<3x!acr.species>, tensor<3xf32>, f32) -> tensor<3x!acr.species>

  acr.return %species_vec_out : tensor<3x!acr.species>
}

acr.reaction @reaction_2(%subvector: tensor<3x!acr.species>, %rate_constant: f32) -> tensor<3x!acr.species> {
  // Define the stoichiometry: CO + O2 -> CO2
  %stoich_coeffs = arith.constant dense<[-1.0, -1.0, 1.0]> : tensor<3xf32>

  // React with the provided rate constant
  %subvector_out = acr.react(%subvector, %stoich_coeffs, %rate_constant) : (tensor<3x!acr.species>, tensor<3xf32>, f32) -> tensor<3x!acr.species>

  acr.return %subvector_out : tensor<3x!acr.species>
}

// ------------------------
// Define the CRN operation using acr.crn_generic with tensor form

acr.crn @reaction_network(%species_vec: tensor<5x!acr.species>, %rate_vec: tensor<2xf32>) -> tensor<5x!acr.species> {
  // For these 5 species: H2, H2O, O2, CO, CO2
  // Stoichiometric matrix encodes both reactions:
  //   2H2 + O2 -> 2H2O  and  CO + O2 -> CO2
  %species_vec_out = acr.crn_generic
      {stoich_matrix = dense<[[-2, -1, 2, 0, 0], [0, -1, 0, -1, 1]]> : tensor<2x5xi32>}
      (%species_vec, %rate_vec)
          : (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>

  acr.return %species_vec_out : tensor<5x!acr.species>
}
```

### Chemical Reaction Network Representation and Exporting to Emulator

The internals of an `acr.crn` can be further canonicalized with an operation that more compactly represents a chemical reaction network `acr.crn_generic`. This operation takes the same inputs, but instead of defining reaction instances in its body, the operation uses a stoichiometric matrix in its attribute to represent the entire reaction network in a compact and efficient manner. The operation can then be exported to a JSON file for further analysis or integration with other tools.

```mlir
%species_out = acr.crn_generic
    {stoich_matrix = dense<[
        // With the following order of species:
        // H2, H2O, O2, CO, CO2
        [-2,-1, 2, 0, 0], // 2H2 + O2 -> 2H2O
        [ 0,-1, 0,-1, 1]  // CO + O2 -> CO2
    ]> : tensor<2x5xi32>}
    (%species_vec, %rate_vec)
        : (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>
```

**Important**: Traditional CRN representations consider reactions as columns, and species as rows in the stoichiometric matrix. In ChemComp we take the transposed approach. In ChemComp `stoich_matrix` attribute has **reactions as rows**, with different **species in each column**. This representation allows for efficient composition of the stoihiometric matrix with individual species vector of `acr.reaction_instance` operations during canonicalization.


### Compilation Flow of Mathemicatical Equations into Chemical Reactions

The process of representing and transforming chemical reactions in ChemComp involves several steps, each corresponding to different levels of abstraction and computation.

#### Mathematical Pattern to `acr.react` / `acr.reaction`

At the most granular level, mathematical patterns are replaced by one or more chemical reactions, represented using the `acr.react` and `acr.reaction` operations. These operations encapsulate the stoichiometry and kinetics of a single reaction.


#### `f(x)` to Multiple `acr.reaction` in an `acr.crn`

A mathematical function `f(x)` can be represented by more complex chemical processes, which involve multiple reactions described using the `acr.crn` operation. This operation encapsulates a network of reactions, allowing for the modeling of interdependent reactions.


#### `acr.crn` to `acr.crn_generic`

The `acr.crn` operation can be transformed or canonicalized into a more generic form, using `acr.crn_generic` operation(s). This operation uses a stoichiometric matrix to represent the entire reaction network in a compact and efficient manner.

The `acr.crn_generic` operation has an attribute that sets up the stoichiometric matrix as a tensor, where rows represent reactions and columns represent species involved in the reactions. It takes as input a tensor of species (concentrations) and a tensor of rates.

```mlir
acr.crn @generic_network(%species_vec: tensor<5x!acr.species>, %rate_vec: tensor<2xf32>) -> tensor<5x!acr.species> {
  %species_vec_out = acr.crn_generic
      {stoich_matrix = dense<[
          // With the following order of species:
          // H2, H2O, O2, CO, CO2
          [-2,-1, 2, 0, 0], // 2H2 + O2 -> 2H2O
          [ 0,-1, 0,-1, 1]  // CO + O2 -> CO2
      ]> : tensor<2x5xi32>}
      (%species_vec, %rate_vec)
          : (tensor<5x!acr.species>, tensor<2xf32>) -> tensor<5x!acr.species>
  acr.return %species_vec_out : tensor<5x!acr.species>
}
```

#### `acr.generic_crn` to JSON File

Finally, the `acr.generic_crn` operation can be exported to a JSON file for further analysis or integration with other tools such as a CRN emulator. The JSON file contains the stoichiometric matrix, species concentrations, and reaction rates, providing a comprehensive representation of the reaction network.

```json
// Chemical Reactions:
// A + B -> C
// C -> D
{
    // Reactions are written in a list.
    "reactions": [
        {
            // // A + B -> C
            "species"   : ["A", "B", "C"],
            "reactants" : [1, 1, 0],
            "products"  : [0, 0, 1],
            "rate": 0.1,
            "name": "R1"
        },
        {
            // C -> D
            "species"   : ["C", "D"],
            "reactants" : [1, 0],
            "products"  : [0, 1],
            "rate": 0.05,
            "name": "R2"
        }
    ],

    // Initial populations
    "populations": {
        "A": 1,
        "B": 1,
        "C": 0,
        "D": 0
    },

    // List of species to track on the simulator, 
    // this is the concentration that should implement our ODE
    "tracked_species": ["D"]
}
```

This flow ensures that chemical reactions are represented accurately and efficiently at different levels of abstractions, facilitating both emulation and verification.


### Dialect rationale

The `acr` dialect operations have in their argument list elements that change and in the body elements that remain constant. For example, for a given reaction, `acr.reaction`, the species and the rate constant are arguments, while the stoichiometry is constant. This separation allows for more flexible and efficient modeling and computation.

When composing mathematical expressions with chemical reactions, fine-grained operations are issued in a pattern-matching fashion. For example, the mathematical equation \( f(x) \) is a pattern that can be matched with a chemical reaction where a species represents \( f(x) \).



### Final Considerations

When building MLIR infrastructure for chemical reactions:
- **Symbolic reasoning**: Allow for symbolic representation of species and reactions to facilitate both simulation and formal verification.
- **Optimization**: Consider how to optimize the representation of large reaction networks, especially for computational chemistry simulations or dynamic reaction systems.
- **Parallelism**: Explore how to parallelize the evaluation of large reaction networks, since many reactions might occur independently.

The infrastructure could potentially integrate with external libraries for chemical informatics (e.g., RDKit, Open Babel) to handle the representation of chemical structures and reactions at a more detailed level.

