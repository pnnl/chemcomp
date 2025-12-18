// RUN: standalone-opt -split-input-file -verify-diagnostics %s | FileCheck %s

// CHECK: func private @species_type(!acr.species) -> ()
func.func private @species_type(!acr.species<>) -> ()

//----

// CHECK: func private @species_type(!acr.species) -> ()
func.func private @vector_species_type(vector<2x!acr.species>) -> ()
