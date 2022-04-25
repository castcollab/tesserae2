from setuptools import setup, Extension
from Cython.Build import cythonize
import Cython.Compiler.Options

Cython.Compiler.Options.annotate = True

setup(
    ext_modules=cythonize([
            Extension("tesserae.hmm", ["src/tesserae/hmm.pyx"],
                      include_dirs=["src/"],
                      extra_compile_args=["-std=c++14", "-O3", "-march=native"])
        ], language_level="3", annotate=True),
)
