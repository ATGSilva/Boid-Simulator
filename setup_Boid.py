from distutils.core import setup
from Cython.Build import cythonize
import numpy


setup(
	include_dirs=[numpy.get_include()],
    ext_modules = cythonize("boid_lib.pyx", compiler_directives={'language_level' : "3", 'infer_types': True})
)