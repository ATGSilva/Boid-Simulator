from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

ext_modules = [
    Extension(
        "boid_lib_mpi_openmp",
        ["boid_lib_mpi_openmp.pyx"],
#  Uncomment the following line if using linux
        extra_compile_args=['-fopenmp','-v'],
		extra_link_args=['-lgomp', '-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/8/'],
		)
]


setup(name="boid_lib_mpi_openmp",
      include_dirs = [np.get_include()],
      ext_modules=cythonize(ext_modules, annotate=True,
                            compiler_directives={'language_level':3,
                                                 'infer_types': True}))