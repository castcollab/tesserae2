import os
import shlex

from setuptools import setup, Extension
from Cython.Build import cythonize
import Cython.Compiler.Options

Cython.Compiler.Options.annotate = True

COMPILE_ARGS = ["-std=c++14", "-O3"]
ARCHITECTURE = os.environ.get('COMPILATION_ARCH')
if ARCHITECTURE:
    arch = shlex.quote(ARCHITECTURE)
    COMPILE_ARGS.append(f"-march={arch}")

setup(
    ext_modules=cythonize([
        Extension("tesserae.hmm", ["src/tesserae/hmm.pyx"],
                  include_dirs=["src/"],
                  extra_compile_args=COMPILE_ARGS)
    ], language_level="3", annotate=True),
)
