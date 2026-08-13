# RUN: %python %s | FileCheck %s

from mlir_standalone.ir import *
from mlir_standalone.dialects import builtin as builtin_d, acr as acr_d

with Context():
    acr_d.register_dialect()
    module = Module.parse(
        """
    %0 = arith.constant 2 : i32
    %2 = acr.stoich_coeff -2 : f32
    """
    )
    # CHECK: %[[C:.*]] = arith.constant 2 : i32
    # CHECK: acr.stoich_coeff -2 : f32
    print(str(module))
