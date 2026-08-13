// RUN: chemcomp-opt -split-input-file -verify-diagnostics %s | FileCheck %s

// Verify that !acr.species parses and prints without angle brackets.

// CHECK: func.func private @species_func(!acr.species)
func.func private @species_func(!acr.species) -> ()

// -----

// CHECK: func.func private @multi_species(!acr.species, !acr.species) -> !acr.species
func.func private @multi_species(%arg0: !acr.species, %arg1: !acr.species) -> !acr.species

// -----

// Verify that tensor<Nx!acr.species> parses and roundtrips.

// CHECK: func.func private @species_tensor(tensor<5x!acr.species>) -> tensor<5x!acr.species>
func.func private @species_tensor(tensor<5x!acr.species>) -> tensor<5x!acr.species>
