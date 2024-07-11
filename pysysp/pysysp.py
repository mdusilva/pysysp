import glob
import os
import warnings
import pprint

import numpy as np
import astropy.io.fits as pyfits
from scipy import interpolate

from .extinction import extlaws

c = 2.99792458e18 #speed of light in Angstron/sec

abzero = -48.60    # magnitude zero points
stzero = -21.10

vega_file = os.path.join(os.path.dirname(__file__), r'alpha_lyr_stis_006.fits')

def Vega(file):
    """Load Vega spectrum"""
    hdulist = pyfits.open(file)
    vega_wavelength = hdulist[1].data.field('Wavelength')
    vega_flux = hdulist[1].data.field('Flux')
    hdulist.close()
    return vega_wavelength, vega_flux

vegawavelength, vegaflux = Vega(vega_file)

filterdir = os.path.join(os.path.dirname(__file__), 'filters')
listfilters = {}
showfilters = {}

def _listfiles(dirname):
    """Create dictionary containing paths to available filters"""
    for root, dirs, files in os.walk(dirname):
        fkey = os.path.basename(root)
        f = []
        for name in files:
            key = os.path.splitext(name)[0]
            listfilters[key] = os.path.join(root, name)
            f.append(key)
        if f:
            showfilters[fkey] = f
        
_listfiles(filterdir)

def add_filter(file_path, name=None, system=None):
    """
    Add a new filter to Pysysp package's library

    Parameters
    ----------
    file_path: string 
        indicates the name of a file with the filter to add.
        It accepts files in plain ASCII 
        with two columns: wavelength and response.
    name: string, default is None
        New filter name (must be different from current names 
        in the current library)
    system: system, default is None
        Filer system, if None it will be saved under "Unknown_system"
    """
    data = np.loadtxt(file_path)
    if len(data.shape) != 2 or data.shape[1] != 2:
        print("File must contain two columns (wavelength, filter response)")
        return
    default_name, fileExtension = os.path.splitext(os.path.basename(file_path))
    if name is None:
        name = default_name
    save_path = filterdir
    if system is None:
        system = "Unknown_system"
    save_path = os.path.join(filterdir, system)
    if not os.path.isdir(save_path):
        os.mkdir(save_path)
    complete_path = os.path.join(save_path, name + ".dat")
    if not os.path.exists(complete_path):
        np.savetxt(complete_path, data)
        print("Filter %s added to filter library" % name)
        _listfiles(filterdir)
    else:
        print("Filter %s could not be added to library because a filter with this name already exists." % name)

#Initialize list of available extincion laws - only Cardelli is implemented at the moment

listlaws = list(extlaws.keys())

def _file_exists(name):
    """Check if a file exists and is accessible"""
    try:
        f = open(name)
        f.close()
        return True
    except IOError:
        return False

class GeneralSpectrum(object):
    """General spectrum class"""

    def __init__(self):
        self.wavelength = None
        self.flux = None

    def loadfits(self, filename=None, wavename='Wavelength', fluxname='Flux'):
        """Load a FITS file containing a spectrum"""
        hdulist = pyfits.open(filename)
        self.wavelength = hdulist[1].data.field(wavename)
        self.flux = hdulist[1].data.field(fluxname)
        hdulist.close()
        
    def loadascii(self, filename=None):
        """Load an ASCII file containing one or two columns"""
        data = np.loadtxt(filename)
        if len(data.shape) == 1:
            self.flux = data
        elif len(data.shape) == 2:
            self.wavelength = data[:,0]
            self.flux = data[:,1]

    def setflux(self, flux):
        self.flux = flux
   
    def setwavelength(self, w):
        self.wavelength = w

class StarSpectrum(GeneralSpectrum):
    """
    Load and operate on spectra
    
    Parameters
    ----------
    file: string 
        indicates the name of a file with the spectrum to load.
        It accepts files with the FITS format or plain ASCII 
        with one (flux) or two columns.
        Default is None

    """
    
    def __init__(self, file=None):
        self.wavelength = None
        self.flux = None
        self._file = file
        if self._file:
            if _file_exists(self._file):
                fname, fextension = os.path.splitext(self._file)
                if fextension == ".fits" or fextension == ".FITS":
                    self.loadfits(self._file)
                else:
                    self.loadascii(self._file)
            else:
                warnings.warn("Warning: Could not find file %s - no spectrum loaded" % (self._file))
                        
    def reflux(self, theta=None):
        """Scale the flux by theta**2"""
        self.flux = self.flux * theta  ** 2.
        
    
    def apmag(self, band, mag='Vega', mzero=0.):
        """
        Compute an apparent magnitude in a given system

        Parameters
        ----------
        band: function or callable 
            a band pass given as a function of wavelength
        mag: {'Vega','AB','ST'}
            name of the magnitude system to use
            Default is 'Vega'
        mstandard: number
            magnitude of the standard star corresponding to
            the zero point (usually 0.) in Vega system
            Default is 0.
        """
        r = self._waverange(band.wavelength, self.wavelength)
        wr = self.wavelength[r]
        fr = self.flux[r]
        if mag == 'AB':
            f = np.trapz(fr * band(wr) * wr, x=wr) / np.trapz(band(wr) * c / wr, x=wr)
            return -2.5 * np.log10(f) + abzero
        elif mag == 'ST':
            f = np.trapz(fr * band(wr) * wr, x=wr) / np.trapz(band(wr) * wr, x=wr)
            return -2.5 * np.log10(f) + stzero
        elif mag == 'Vega':
            vegar = self._waverange(band.wavelength, vegawavelength)
            vegawr = vegawavelength[vegar]
            vegafr = vegaflux[vegar]
            f = np.trapz(fr * band(wr) * wr, x=wr)
            vega_f = np.trapz(vegafr * band(vegawr) * vegawr, x=vegawr)
            return -2.5 * np.log10(f) + 2.5 * np.log10(vega_f) + mzero
        else:
            raise ValueError('Magnitude system is not a valid choice, check input string')
        
    def extinction(self, band, law='cardelli', **kwargs):
        """
        Computes extinction in band assuming an extinction law.
        Input the desired band and law. The law should be a callable as 
        function of wavelength.
        
        Aditional keyword arguments may be passed to the extinction law.
        """
        r = self._waverange(band.wavelength, self.wavelength)
        wr = self.wavelength[r]
        fr = self.flux[r]
        if law in extlaws:
            f = np.trapz(fr * band(wr) * wr * 10.**(-0.4 * extlaws[law](wr,**kwargs)), \
            x=wr) / np.trapz(fr * band(wr) * wr, x=wr)
        else:
            raise ValueError('Extinction law is not a valid choice, check input string')
        return -2.5 * np.log10(f)
            
    def _waverange(self, wb, ws):
        wbmax = np.max(wb)
        wbmin = np.min(wb)
        wsmax = np.max(ws)
        wsmin = np.min(ws)
        if wbmin < wsmin or wbmax > wsmax:
            raise ValueError('Bandpass and spectrum must completely overlap.')
        else:
            range = np.where(np.logical_and(ws <= wbmax, ws >= wbmin))
            return range
            
class BandPass(GeneralSpectrum):
    """
    Bandpasses photonic response curves
    
    Parameters
    ----------
    band: string 
        indicates the name of a filter in the internal list (see showfilters function)
          or a filename of a file containing two columns, one with the wavelength and 
          one with the normalized response
    smt: string 
        indicates the type of smoothing to perform on the repsonse curve. Uses same
        names as Scipy's interp1d ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')
        Default is slinear
    
    """
    
    def __init__(self, band=None, smt='slinear'):
        """Initialize a bandpass"""
        self.wavelength = None
        self.flux = None
        self.name = None
        self.interpolated_response = None
        self._smt = smt
        #Check if input is a band part of our list or an existing file
        if band:
            if band in listfilters:
                self.loadascii(listfilters[band])
                if self._smt in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
                    self.smooth(kind=self._smt)
                else:
                    raise ValueError('Unknown type of smoothing')

            elif _file_exists(band):
                self.loadascii(band)
                if self._smt in ['linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic']:
                    self.smooth(kind=self._smt)
                else:
                    raise ValueError('Unknown type of smoothing')

            else:
                warnings.warn("Warning: Could not find filter or file %s " % (band))
        
    def loadascii(self, filename=None):
        data = np.loadtxt(filename)
        self.name, fileExtension = os.path.splitext(os.path.basename(filename))
        self.wavelength = data[:,0]
        self.response = data[:,1]

    def smooth(self, kind='linear'):
        if interpolate is not None:
            self.interpolated_response = interpolate.interp1d(self.wavelength, self.response, kind=kind, bounds_error=False, fill_value=0.)
        else:
            warnings.warn("Warning: Smoothing requires Scipy, returning a linear interpolation instead.")
            self.interpolated_response = lambda w: np.interp(w, self.wavelength, self.response, left=0., right=0.)
            
    def __call__(self, w):
        return self.interpolated_response(w)
