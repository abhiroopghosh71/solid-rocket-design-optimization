from setuptools import setup
from setuptools.extension import Extension

import numpy as np
from Cython.Build import cythonize

approcket_extension = Extension(
    name="approcket",
    sources=["approcket.pyx"],
    libraries=["approcket"],
    library_dirs=["librocket"],
    include_dirs=["librocket", np.get_include()]
)
setup(
    name="approcket",
    ext_modules=cythonize([approcket_extension], gdb_debug=True, compiler_directives={'language_level': "3"})
)
