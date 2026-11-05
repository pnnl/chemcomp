# -*- Python -*-

import os
import platform
import re
import subprocess
import sys
import tempfile

import lit.formats
import lit.util

from lit.llvm import llvm_config
from lit.llvm.subst import ToolSubst
from lit.llvm.subst import FindTool

# Configuration file for the 'lit' test runner.

# name: The name of this test suite.
config.name = "CHEMCOMP"

config.test_format = lit.formats.ShTest(not llvm_config.use_lit_shell)

# suffixes: A list of file extensions to treat as test files.
config.suffixes = [".mlir"]

# test_source_root: The root path where tests are located.
config.test_source_root = os.path.dirname(__file__)

# test_exec_root: The root path where tests should be run.
config.test_exec_root = os.path.join(config.chemcomp_obj_root, "test")

config.substitutions.append(("%PATH%", config.environment["PATH"]))
config.substitutions.append(("%shlibext", config.llvm_shlib_ext))

llvm_config.with_system_environment(["HOME", "INCLUDE", "LIB", "TMP", "TEMP"])

llvm_config.use_default_substitutions()

# excludes: A list of directories to exclude from the testsuite. The 'Inputs'
# subdirectories contain auxiliary inputs for various tests in their parent
# directories.
config.excludes = ["Inputs", "Examples", "CMakeLists.txt", "README.txt", "LICENSE.txt",
                   "ACRN", "Transform"]

# The Base Modelica (MARCO) frontend is only built with CHEMCOMP_ENABLE_MARCO.
if config.enable_marco:
    config.available_features.add("marco")

# test_exec_root: The root path where tests should be run.
config.test_exec_root = os.path.join(config.chemcomp_obj_root, "test")
config.chemcomp_tools_dir = os.path.join(config.chemcomp_obj_root, "bin")
config.chemcomp_libs_dir = os.path.join(config.chemcomp_obj_root, "lib")

config.substitutions.append(("%chemcomp_libs", config.chemcomp_libs_dir))

# Tweak the PATH to include the tools dir.
llvm_config.with_environment("PATH", config.llvm_tools_dir, append_path=True)

tool_dirs = [config.chemcomp_tools_dir, config.llvm_tools_dir]
tools = [
    "chemcomp-tblgen",
    "mlir-opt",
    "chemcomp-capi-test",
    "chemcomp-opt",
    "chemcomp-translate",
]

llvm_config.add_tool_substitutions(tools, tool_dirs)

llvm_config.with_environment(
    "PYTHONPATH",
    [
        os.path.join(config.mlir_obj_dir, "python_packages", "chemcomp"),
    ],
    append_path=True,
)

# The ODE-to-network examples check their networks numerically with the
# reservoir emulator when its Python package is importable (the default pixi
# environment). %chemcomp_python runs the interpreter lit itself runs under,
# with -s: lit does not pass PYTHONNOUSERSITE on, and packages in
# ~/.local must not shadow the environment's.
config.substitutions.append(("%chemcomp_python", sys.executable + " -s"))
config.substitutions.append(
    (
        "%check_realization",
        os.path.join(
            os.path.dirname(config.test_source_root),
            "tools",
            "chemcomp-reservoir-emulator",
            "scripts",
            "check_realization.py",
        ),
    )
)
if (
    subprocess.run(
        [sys.executable, "-s", "-c", "import chemical_reservoir, scipy"],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode
    == 0
):
    config.available_features.add("emulator")
