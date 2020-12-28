from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy

ext_modules = [
    Extension(
        "boid_lib_mpi",
        ["boid_lib_mpi.pyx"],
#  Uncomment the following line if using linux
        extra_compile_args=['-fopenmp','-v'],
		extra_link_args=['-lgomp', '-Wl,-rpath,/usr/local/opt/gcc/lib/gcc/8/'],
		)
]

setup(
	include_dirs=[numpy.get_include()],
    ext_modules = cythonize(ext_modules, annotate=True,
                            compiler_directives={'language_level':3,
                                                 'infer_types': True})
)