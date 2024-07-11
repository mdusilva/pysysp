"""
pysysp: Synthetic Stellar Photometry for Python

This is a set of tools to compute synthetic photometry in a simple way,
ideal to integrate in larger projects.
The inputs are photonic response functions for the desired phtometric bands and 
stellar spectra. The units used in this package are: erg cm^-2 s^-1 Ang^-1 for fluxes; 
Ang^-1 for wavelength. Note that this implies units of cm Ang^-1 for the speed of light.
These units must be used in all inputs for the package to operate properly.
It is possible to compute magnitudes in Vega and AB systems, and interstellar 
extinction, for any given photometric passband, given an extinction law, which must 
be a callable object function of wavelength only. The only extinction law provided 
is the one by Cardelli, Clayton & Mathis (1989), but others can be easily added if 
needed.
At the moment only the bandpasses for Bessel's UBVRI and Gaia's G filters are provided.

References: Casagrande & VandenBerg (2014); Bessel & Murphy (2012)
Acknowledgments:  Vega CALSPEC spectrum alpha_lyr_stis_006
                http://www.stsci.edu/hst/observatory/cdbs/calspec.html
                  photonic response functions by Bessel & Murphy (2012)
                  for the U, B, V, R, I filters




-------------------------
Usage example:

>>> import pysysp

#load a spectrum from a fits file
>>> vega = pysysp.StarSpectrum('alpha_lyr_stis_006.fits')

#print the flux
>>> vega.flux
array([  1.23810533e-17,   1.67559561e-17,   1.78002369e-17, ...,
         1.40140733e-19,   1.38734358e-19,   1.26490663e-19], dtype=float32)

#Check available filters
>>> pysysp.showfilters
{'Gaia': ['G'], 'Bessel': ['B', 'I', 'R', 'U', 'V']}
         
#load the photonic response for the bands B and V with linear smoothing
>>> V=pysysp.BandPass('V',smt='linear')
>>> B=pysysp.BandPass('B',smt='linear')

#You can also load any ascii file with two columns (first wavelength, second filter response)
>>> B2 = pysysp.BandPass('./filters/B.dat')

#You can also permanetly add a new filter to the libary of available filters
>>> pysysp.add_filter('test.txt', name='T')

#compute the V magnitude in Vega system
#The magnitude of Vega in V band is not 0, so in order to
#compute the correct value we should correct the zero point
>>> vega.apmag(V,mag='Vega',mzero=0.03)
0.029999999999999999

#Same in AB system
>>> vega.apmag(V,mag='AB')
0.0062420242892642364


#Compute E(B-V) for a (monochromatic) extinction Av=0.5
>>> ab=vega.extinction(B,law='cardelli',A=0.5,Rv=3.1)
>>> av=vega.extinction(V,law='cardelli',A=0.5,Rv=3.1)
>>> ab-av
0.15954421668471008

"""

from .pysysp import StarSpectrum, BandPass, showfilters, listlaws, add_filter

__version__ = '2.0.1'