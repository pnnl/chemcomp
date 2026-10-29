// The MARCO Base Modelica dialect is registered in chemcomp-opt (built with
// CHEMCOMP_ENABLE_MARCO=ON), so ODE models in the form produced by the MARCO
// frontend parse, verify and print back unchanged.

// RUN: chemcomp-opt %s | chemcomp-opt | FileCheck %s

// Exponential decay, der(x) = -k*x, as emitted by the MARCO frontend: the
// derivative is still a bmodelica.der operation.

// CHECK-LABEL: bmodelica.model @Decay {
// CHECK:         bmodelica.variable @x : !bmodelica.variable<!bmodelica.real>
// CHECK:         bmodelica.variable @k : !bmodelica.variable<!bmodelica.real, parameter>
// CHECK:         bmodelica.binding_equation @k {
// CHECK:           bmodelica.constant #bmodelica<real 5.000000e-01>
// CHECK:         %[[T:.*]] = bmodelica.equation_template inductions = [] {
// CHECK:           %[[X:.*]] = bmodelica.variable.get @x : !bmodelica.real
// CHECK:           %[[DX:.*]] = bmodelica.der %[[X]] : !bmodelica.real -> !bmodelica.real
// CHECK:           %[[K:.*]] = bmodelica.variable.get @k : !bmodelica.real
// CHECK:           %[[KX:.*]] = bmodelica.mul %[[K]], %[[X]]
// CHECK:           %[[RHS:.*]] = bmodelica.neg %[[KX]]
// CHECK:           bmodelica.equation_sides
// CHECK:         bmodelica.dynamic {
// CHECK-NEXT:      bmodelica.equation_instance %[[T]]
bmodelica.model @Decay {
    bmodelica.variable @x : !bmodelica.variable<!bmodelica.real>
    bmodelica.variable @k : !bmodelica.variable<!bmodelica.real, parameter>

    bmodelica.binding_equation @k {
        %0 = bmodelica.constant #bmodelica<real 0.5>
        bmodelica.yield %0 : !bmodelica.real
    }

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !bmodelica.real
        %dx = bmodelica.der %x : !bmodelica.real -> !bmodelica.real
        %k = bmodelica.variable.get @k : !bmodelica.real
        %kx = bmodelica.mul %k, %x : !bmodelica.real, !bmodelica.real -> !bmodelica.real
        %rhs = bmodelica.neg %kx : !bmodelica.real -> !bmodelica.real
        %l = bmodelica.equation_side %dx : tuple<!bmodelica.real>
        %r = bmodelica.equation_side %rhs : tuple<!bmodelica.real>
        bmodelica.equation_sides %l, %r : tuple<!bmodelica.real>, tuple<!bmodelica.real>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}

// After derivatives-materialization the model carries a derivatives map and
// der_x is an ordinary variable.

// CHECK-LABEL: bmodelica.model @Materialized der = [<@x, @der_x>] {
// CHECK:         bmodelica.variable @der_x : !bmodelica.variable<!bmodelica.real>
// CHECK:         bmodelica.variable.get @der_x
bmodelica.model @Materialized der = [<@x, @der_x>] {
    bmodelica.variable @x : !bmodelica.variable<!bmodelica.real>
    bmodelica.variable @der_x : !bmodelica.variable<!bmodelica.real>

    %t0 = bmodelica.equation_template inductions = [] {
        %x = bmodelica.variable.get @x : !bmodelica.real
        %dx = bmodelica.variable.get @der_x : !bmodelica.real
        %rhs = bmodelica.neg %x : !bmodelica.real -> !bmodelica.real
        %l = bmodelica.equation_side %dx : tuple<!bmodelica.real>
        %r = bmodelica.equation_side %rhs : tuple<!bmodelica.real>
        bmodelica.equation_sides %l, %r : tuple<!bmodelica.real>, tuple<!bmodelica.real>
    }

    bmodelica.dynamic {
        bmodelica.equation_instance %t0
    }
}
