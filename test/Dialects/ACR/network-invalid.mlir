// RUN: chemcomp-opt %s -split-input-file -verify-diagnostics

acr.network @N {
  acr.species @x
  // expected-error @below {{'@y' is not an acr.species of this network}}
  acr.rxn @x -> @y
}

// -----

acr.network @N {
  acr.species @x
  // expected-error @below {{reactants and products are the same complex}}
  acr.rxn @x -> @x
}

// -----

acr.network @N {
  acr.species @x
  // expected-error @below {{rate constants must not be negative}}
  acr.rxn @x -> none rate -1.0
}

// -----

acr.network @N {
  acr.species @x
  // expected-error @below {{species listed twice in a complex}}
  acr.rxn @x + @x -> none
}

// -----

acr.network @N {
  acr.species @x
  // expected-error @below {{stoichiometric counts must be positive}}
  acr.rxn 0 @x -> none
}
