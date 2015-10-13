from distutils.core import setup
from distutils.core import Extension
from Cython.Build import cythonize


setup ( ext_modules = cythonize( [Extension("cphotometry", ["cphotometry.pyx"], libraries=["gsl","cblas"])] ) )

