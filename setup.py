from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

examples_extension = Extension(
    name="linear_interp",
    sources=["linear_interp.pyx", "lib/bin_interp.c"],
    # sources=["linear_interp.pyx"],
    language="c",
    # libraries=["bin_interp"],
    # libraries=["src"],
    include_dirs=["lib", numpy.get_include()],
    library_dirs=["lib"],
)

setup(
    name="linear_interp",
    ext_modules=cythonize([examples_extension])
)