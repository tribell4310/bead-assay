"""
build.py
"""

import distutils.core
import Cython.Build
import numpy

distutils.core.setup(ext_modules = Cython.Build.cythonize("bead_mask_analyze_cython_modules.pyx"), include_dirs=[numpy.get_include()])