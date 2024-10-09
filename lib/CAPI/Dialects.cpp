//===- Dialects.cpp - CAPI for dialects -----------------------------------===//
//
// Part of the LLVM Project, under the Apache License v2.0 with LLVM Exceptions.
// See https://llvm.org/LICENSE.txt for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//===----------------------------------------------------------------------===//

#include "Standalone-c/Dialects.h"

#include "Standalone/Dialect/ACR/ACRDialect.h"
#include "Standalone/StandaloneDialect.h"
#include "mlir/CAPI/Registration.h"

MLIR_DEFINE_CAPI_DIALECT_REGISTRATION(Standalone, standalone,
                                      mlir::standalone::StandaloneDialect)

MLIR_DEFINE_CAPI_DIALECT_REGISTRATION(ACR, acr, mlir::acr::ACRDialect)