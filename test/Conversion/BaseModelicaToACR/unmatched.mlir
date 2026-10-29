// Models outside the supported subset, or matching no known network, are
// diagnosed rather than silently dropped.

// RUN: chemcomp-opt %s --convert-bmodelica-to-acr -split-input-file -verify-diagnostics

!r = !bmodelica.real

// Dimerization-like decay der(x) = -k*x^2 is polynomial, but not in the table.
// expected-error @below {{no known chemical reaction network matches the ODE system of model 'Quadratic'}}
// expected-note @below {{ODE system: der(x) = -k*x^2}}
bmodelica.model @Quadratic {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @k : !bmodelica.variable<!r, parameter>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        %k = bmodelica.variable.get @k : !r
        %0 = bmodelica.mul %x, %x : !r, !r -> !r
        %1 = bmodelica.mul %k, %0 : !r, !r -> !r
        %2 = bmodelica.neg %1 : !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// -----

!r = !bmodelica.real

// Lotka-Volterra shape, but predation does not conserve mass: the two
// interaction terms have different rate constants, so no X + Y -> 2Y.
// expected-error @below {{no known chemical reaction network matches the ODE system of model 'Unbalanced'}}
// expected-note @below {{ODE system: der(x) = -b1*x*y + x; der(y) = b2*x*y - y}}
bmodelica.model @Unbalanced {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @y : !bmodelica.variable<!r>
    bmodelica.variable @b1 : !bmodelica.variable<!r, parameter>
    bmodelica.variable @b2 : !bmodelica.variable<!r, parameter>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %dx = bmodelica.der %x : !r -> !r
        %b = bmodelica.variable.get @b1 : !r
        %0 = bmodelica.mul %b, %x : !r, !r -> !r
        %1 = bmodelica.mul %0, %y : !r, !r -> !r
        %2 = bmodelica.sub %x, %1 : !r, !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    %t1 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %dy = bmodelica.der %y : !r -> !r
        %b = bmodelica.variable.get @b2 : !r
        %0 = bmodelica.mul %b, %x : !r, !r -> !r
        %1 = bmodelica.mul %0, %y : !r, !r -> !r
        %2 = bmodelica.sub %1, %y : !r, !r -> !r
        %l = bmodelica.equation_side %dy : tuple<!r>
        %r = bmodelica.equation_side %2 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
        bmodelica.equation_instance %t1
    }
}

// -----

!r = !bmodelica.real

// Non-polynomial right-hand sides are outside the supported subset.
bmodelica.model @Trigonometric {
    bmodelica.variable @x : !bmodelica.variable<!r>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %dx = bmodelica.der %x : !r -> !r
        // expected-error @below {{unsupported operation 'bmodelica.sin' in an ODE expression}}
        %0 = bmodelica.sin %x : !r -> !r
        %l = bmodelica.equation_side %dx : tuple<!r>
        %r = bmodelica.equation_side %0 : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// -----

!r = !bmodelica.real

// Algebraic equations (no derivative) are not ODEs.
bmodelica.model @Algebraic {
    bmodelica.variable @x : !bmodelica.variable<!r>
    bmodelica.variable @y : !bmodelica.variable<!r>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !r
        %y = bmodelica.variable.get @y : !r
        %l = bmodelica.equation_side %y : tuple<!r>
        %r = bmodelica.equation_side %x : tuple<!r>
        bmodelica.equation_sides %l, %r : tuple<!r>, tuple<!r>
    }

    bmodelica.dynamic {
        // expected-error @below {{algebraic equations are not supported; every equation must define the derivative of a state}}
        bmodelica.equation_instance %t0
    }
}
