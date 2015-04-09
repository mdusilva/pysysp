import os
import re

from setuptools import setup

vre = re.compile("__version__ = \'(.*?)\'")
m = open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "pysysp", "__init__.py")).read()
version = vre.findall(m)[0]

setup(
      name = "pysysp",
      version = version,
      packages = ["pysysp"],
      description = "A package for Synthetic Stellar Photometry",
      long_description=open("README.rst").read() + "\n\n"
                + "Changelog\n"
                + "---------\n\n"
                + open("HISTORY.rst").read(),

      author = "Manuel Silva",
      author_email = "madusilva@gmail.com",
      url='https://github.com/mdusilva/pysysp'
      license="GPLv2",
      classifiers=[
          "Intended Audience :: Developers",
          "Intended Audience :: Science/Research",
          "License :: OSI Approved :: GNU General Public License (GPL)",
          "Operating System :: OS Independent",
          "Programming Language :: Python",
          "Topic :: Scientific/Engineering :: Astronomy"

      ],
      install_requires = ["numpy", "pyfits"],
      package_data = {
          '' : ['alpha_lyr_stis_006.fits', 'filters/Bessel/*.dat','filters/Gaia/*.dat']
      },
      zip_safe=False
)