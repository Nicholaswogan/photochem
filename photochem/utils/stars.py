import numpy as np
import requests
import tempfile
from astropy.io import fits
import numba as nb
from scipy import constants as const
import yaml
import h5py
from photochem_clima_data import DATA_DIR
from astroquery.mast import Observations
from urllib.parse import quote

# Relative imports
from .youngsun import youngsun
from .._clima import rebin, rebin_with_errors

###
### Some utilities for dealing with spectra
###

@nb.njit()
def stefan_boltzmann(T):
    """The Stefan-Boltzmann law

    Parameters
    ----------
    T : float
        Temperature (K)

    Returns
    -------
    float
        Energy flux (W/m^2)
    """    
    return const.sigma*T**4

@nb.njit()
def stefan_boltzmann_inverse(F):
    """Inverse of the Stefan-Boltzmann law

    Parameters
    ----------
    F : float
        Energy flux (W/m^2)

    Returns
    -------
    float
        Temperature (K)
    """    
    return (F/const.sigma)**(1/4)

@nb.njit()
def equilibrium_temperature(stellar_radiation, bond_albedo):
    """Equilibrium temperature of a planet.

    Parameters
    ----------
    stellar_radiation : float
        Incident energy at a planet in W/m^2
    bond_albedo : float
        Bond albedo

    Returns
    -------
    float
        The equilibrium temperature (K)
    """    
    T_eq = ((stellar_radiation*(1.0 - bond_albedo))/(4.0*const.sigma))**(0.25)
    return T_eq 

@nb.njit()
def equilibrium_temperature_inverse(Teq, bond_albedo):
    """Inverse of the equation for equilibrium temperature.

    Parameters
    ----------
    Teq : float
        The equilibrium temperature (K)
    bond_albedo : float
        Bond albedo

    Returns
    -------
    float
        Incident energy at a planet in W/m^2
    """    
    stellar_radiation = 4.0*const.sigma*Teq**4/(1.0 - bond_albedo)
    return stellar_radiation

@nb.njit()
def blackbody_cgs(T, lam):
    """Blackbody flux in cgs units in per unit wavelength (erg/cm^2/s/cm/sr)

    Parameters
    ----------
    T : float
        Temperature (K)
    lam : ndarray[ndim=1,double]
        Wavelength (cm)
    
    Returns
    -------
    ndarray[ndim=1,double]
        The blackbody flux at the input wavelengths in erg/cm^2/s/cm/sr.
    """
    h = const.h*1e7 # erg s 
    c = const.c*1e2 # cm/s
    k = const.k*1e7 # erg / K
    B = ((2.0*h*c**2.0)/(lam**5.0))*(1.0/(np.exp((h*c)/(lam*k*T)) - 1.0))
    return B

@nb.njit()
def blackbody(T, wv):
    """Blackbody flux

    Parameters
    ----------
    T : float
        Temperature (K)
    wv : ndarray[ndim=1,double]
        Wavelength (nm)

    Returns
    -------
    ndarray[ndim=1,double]
        The blackbody flux at the input wavelengths in mW/m^2/s/nm/sr.
    """    

    # nm * (m/1e9 nm) * (1e2cm/ m) = cm
    wv_cgs = wv*(1/1e9)*(1e2/1)
    B_cgs = blackbody_cgs(T, wv_cgs)
    # erg/cm^2/s/cm/sr * (1 W/ 1e7 erg) * (1e3 mW/ 1 W) * (1e4 cm^2/ 1 m^2) * 
    # (1e2 cm/ 1 m) * (1 m / 1e9 nm) = mW/m^2/s/nm/sr
    B = B_cgs*(1/1e7)*(1e3/1)*(1e4/1)*(1e2/1)*(1/1e9)
    return B

@nb.njit()
def grid_at_resolution(min_wv, max_wv, R):
    """Computes a grid of bins at a given resolution `R`

    Parameters
    ----------
    min_wv : float
        Minimum bin extent
    max_wv : float
        Maximum bin extent
    R : float
        in resolution (dlam = lam/R)

    Returns
    -------
    ndarray[ndim=1,double]
        Edges of the wavelength grid
    """    
    wavl = [min_wv]
    while wavl[-1] < max_wv:
        dlam = wavl[-1]/R
        wavl.append(wavl[-1]+dlam)
    wavl[-1] = max_wv
    return np.array(wavl)

@nb.njit()
def make_bins(wavs):
    """Given a series of wavelength points, find the edges
    of corresponding wavelength bins.

    Parameters
    ----------
    wavs : ndarray[ndim=1,double]
        Input wavelengths

    Returns
    -------
    ndarray[ndim=1,double]
        Edges of the wavelength grid
    """
    if np.any(wavs[1:] <= wavs[:-1]):
        raise ValueError('Wavelengths `wavs` are not strictly increasing')
    edges = np.zeros(wavs.shape[0]+1)
    edges[0] = wavs[0] - (wavs[1] - wavs[0])/2
    edges[-1] = wavs[-1] + (wavs[-1] - wavs[-2])/2
    edges[1:-1] = (wavs[1:] + wavs[:-1])/2
    return edges

@nb.njit()
def energy_in_spectrum(wv, F):
    """The total energy in a spectrum

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)

    Returns
    -------
    float
        Energy in W/m^2
    """    
    wavl = make_bins(wv)
    return 1e-3*np.sum(F[:]*(wavl[1:]-wavl[:-1])) # W/m^2

def scale_spectrum_to_planet(wv, F, Teq=None, stellar_flux=None):
    """Scales a stellar spectrum so that it has the correct incident
    flux at a planet. You must supply either a `Teq` or a `stellar_flux`

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None

    Returns
    -------
    ndarray[ndim=1,double]
        The stellar spectrum rescaled so that it has the proper total
        bolometric flux at the planet (mW/m^2/nm).
    """    

    if Teq is None and stellar_flux is None:
        raise ValueError('Either `Teq` or `stellar_flux` must be supplied as inputs')
    elif Teq is not None and stellar_flux is not None:
        raise ValueError('Only one of the following inputs should be supplied: `Teq` and `stellar_flux`')
    if np.any(wv[1:] <= wv[:-1]):
        raise ValueError('Wavelengths `wv` are not strictly increasing')
    if wv[0] > 100:
        raise ValueError('Wavelengths `wv` must extend below 100 nm')
    if wv[-1] < 50e3:
        raise ValueError('Wavelengths `wv` must extend above 50 um')
    
    if Teq is not None:
        stellar_flux = equilibrium_temperature_inverse(Teq, 0.0)

    factor = stellar_flux/energy_in_spectrum(wv, F)
    F *= factor
    
    return F

@nb.njit()
def append_blackbody_to_stellar_spectrum(wv, F, Teff, wv_end=100e3, nwb=1000):
    """Appends a blackbody to a stellar spectrum

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)
    Teff : float
        Stellar effective temperature
    wv_end : float, optional
        Ending wavelength for the blackbody (nm), by default 100e3
    nwb : int, optional
        Number of blackbody bins, by default 1000

    Returns
    -------
    wv_new : ndarray[ndim=1,double]
        New wavelengths of stellar spectrum (nm)
    F_new : ndarray[ndim=1,double]
        New flux of the stellar spectrum (mW/m^2/nm)
    """    

    if wv[-1] < wv_end:
        wv1 = np.logspace(np.log10(wv[-1]*(1 + 1e-5)), np.log10(wv_end), nwb)
        F1 = blackbody(Teff, wv1)*np.pi
        factor = F[-1]/(blackbody(Teff, wv[-1])*np.pi)
        F1 *= factor
        wv_new = np.append(wv, wv1)
        F_new = np.append(F, F1)
    
    return wv_new, F_new

def photochem_spectrum_string(wv, F, Teq=None, stellar_flux=None, scale_to_planet=True, fmt='{:20}'):
    """Rescales a stellar spectrum to a planet, then converts it a string
    representing a photochem stellar flux file.

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None
    scale_to_planet : bool, optional
        If True, then the stellar flux will be rescaled to the planet, be default True.
    fmt : str, optional
        Format string, by default '{:20}'

    Returns
    -------
    str
        A photochem stellar flux file as a string.
    """    

    if scale_to_planet:
        F_save = scale_spectrum_to_planet(wv, F, Teq, stellar_flux)
    else:
        F_save = F

    flux_str = ""
    flux_str += fmt.format('Wavelength(nm)')
    flux_str += fmt.format('SolarFlux(mW/m^2/nm)')
    flux_str += '\n'
    for i in range(wv.shape[0]):
        flux_str += fmt.format('%.8e'%wv[i])
        flux_str += fmt.format('%.8e'%F_save[i])
        flux_str += '\n'

    return flux_str

def save_photochem_spectrum(wv, F, outputfile, Teq=None, stellar_flux=None, scale_to_planet=True, fmt='{:20}'):
    """Rescales a stellar spectrum to a planet, if desired, then saves the spectrum in the photochem
    format.

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)
    outputfile : str
        Path to the output stellar file.
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None
    scale_to_planet : bool, optional
        If True, then the stellar flux will be rescaled to the planet, be default True.
    fmt : str, optional
        Format string, by default '{:20}'
    """    

    outstr = photochem_spectrum_string(wv, F, Teq, stellar_flux, scale_to_planet, fmt)
    with open(outputfile,'w') as f:
        f.write(outstr)

def rebin_to_needed_resolution(wv, F):
    """Rebins a stellar spectrum to 4x the resolution needed by the photochemical
    and climate models.

    Parameters
    ----------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)

    Returns
    -------
    wv_new : ndarray[ndim=1,double]
        New wavelengths of stellar spectrum (nm)
    F_new : ndarray[ndim=1,double]
        New flux of the stellar spectrum (mW/m^2/nm)
    """    

    # Get UV bins, and solar radiative transfer bins
    with h5py.File(DATA_DIR+'/xsections/bins.h5','r') as f:
        wavl_uv = f['wavl'][:].astype(np.double)
    with h5py.File(DATA_DIR+'/kdistributions/bins.h5','r') as f:
        wavl_sol = f['sol_wavl'][:].astype(np.double)*1e3 # to nm
    
    # Combine UV and solar bins
    wavl_new = np.sort(np.append(wavl_uv,wavl_sol))
    wavl_new = np.unique(wavl_new)
    # Double the resolution
    wavl_new = np.sort(np.append(wavl_new, (wavl_new[1:] + wavl_new[:-1])/2))
    # Double the resolution again
    wavl_new = np.sort(np.append(wavl_new, (wavl_new[1:] + wavl_new[:-1])/2))

    # Mid-points between bins
    wv_new = (wavl_new[1:] + wavl_new[:-1])/2
    # Add edges of grid to ensure we get that energy
    wv_new = np.append(wavl_new[0], wv_new)
    wv_new = np.append(wv_new, wavl_new[-1])
    wavl_new = make_bins(wv_new)

    # Rebin the spectrum to the needed resolution
    wavl = make_bins(wv)
    F_new = rebin(wavl.copy(), F.copy(), wavl_new.copy())

    return wv_new, F_new

###
### The Sun's spectrum through time.
###

def solar_spectrum(outputfile=None, age=0.0, append_blackbody=True, Teq=None, stellar_flux=None, scale_before_age=True,
                   needed_resolution=True):
    """The Sun's spectrum throughout all time. For the Modern Sun, we use the
    Thuillier et al. (2004) (https://doi.org/10.1029/141GM13), ATLAS 3 reference
    spectrum. Then we scale the spectrum into the past/future with the Claire et al. (2012)
    YoungSun routine.

    Parameters
    ----------
    outputfile : str, optional
        If supplied, this is the output file name, by default None
    age : float, optional
        Age of the sun in billions of years ago, by default 0.0
    append_blackbody : bool, optional
        If True, then a blackbody is appended, by default True
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None
    scale_before_age : bool, optional
        If True, then the spectrum is scaled before the `youngsun` routine, 
        is applied. If False, then scaling occurs afterword, by default True.
    needed_resolution : bool, optional
        If True, then the spectrum is rebinned to a resolution 4x higher than
        What is used by the photochemical and climate models which should be
        an adequately high resolution, by default True.

    Returns
    -------
    wv : ndarray[ndim=1,double]
        Wavelengths of stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum at Earth (mW/m^2/nm)
    """    

    # Download the spectrum
    url = 'https://sbuv.gsfc.nasa.gov/solar/reference_spectra/ATLAS1_2004.txt'
    response = requests.get(url)

    # Error if download didn't work.
    if response.status_code != 200:
        raise Exception("Failed to download the Modern Sun's spectrum")
    
    # Get wavelength and flux
    lines = response.content.decode().split('\n')
    wv = np.array([float(line.split()[0]) for line in lines if len(line) > 0])
    F = np.array([float(line.split()[1]) for line in lines if len(line) > 0])

    # Append a blackbody to long wavelengths, if desired. We use an
    # effective temperature of 5778 K for the Sun.
    if append_blackbody:
        wv, F = append_blackbody_to_stellar_spectrum(wv, F, 5778)

    # Rescale to stellar_flux or Teq
    if scale_before_age and (stellar_flux is not None or Teq is not None):
        if not append_blackbody:
            raise ValueError('You must append a blackbody to rescale the spectrum')
        F = scale_spectrum_to_planet(wv, F, Teq, stellar_flux)

    # Change the age of the spectrum
    if age != 0.0:
        fluxmult = youngsun(age, wv*10)
        F *= fluxmult

    # Rescale to stellar_flux or Teq
    if not scale_before_age and (stellar_flux is not None or Teq is not None):
        if not append_blackbody:
            raise ValueError('You must append a blackbody to rescale the spectrum')
        F = scale_spectrum_to_planet(wv, F, Teq, stellar_flux)

    # Downbin to needed resolution
    if needed_resolution:
        wv, F = rebin_to_needed_resolution(wv, F)

    # Save the file, if desired
    if outputfile is not None:
        save_photochem_spectrum(wv, F, outputfile, scale_to_planet=False)
        
    return wv, F

###
### HAZMAT spectra (https://archive.stsci.edu/hlsp/hazmat)
###

def hazmat_spectrum(star_name, model='model', outputfile=None, Teq=None, stellar_flux=None, needed_resolution=True):
    """Downloads a HAZMAT spectrum (https://archive.stsci.edu/hlsp/hazmat), 
    then rescale the spectrum so that it has a total bolometric insolation at 
    a planet consistent with `Teq` or `stellar_flux`. Finally, the spectrum can 
    be saved to `outputfile` in photochem format.

    Parameters
    ----------
    star_name : str
        Name of the HAZMAT star
    model : str
        Model group (for TRAPPIST-1), either "1a", "2a", or "2b", or "model" for GJ stars.
    outputfile : str, optional
        If not None, then this is the output filename, be default None.
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None
    needed_resolution : bool, optional
        If True, then the spectrum is rebinned to a resolution 4x higher than
        What is used by the photochemical and climate models which should be
        an adequately high resolution, by default True.
    """

    if star_name.lower() == 'trappist-1' and model not in ['1a','2a','2b']:
        raise ValueError('TRAPPIST-1 is only compatible with model 1a, 2a, or 2b')

    base_url = 'http://archive.stsci.edu/hlsps/hazmat/hlsp_hazmat_phoenix_synthspec_'
    url = base_url+star_name.lower()+'_'+model+'_v1_fullres.fits'
    response = requests.get(url)

    if response.status_code != 200:
        raise Exception('Failed to download '+star_name+' from HAZMAT')
    
    with tempfile.TemporaryFile() as f:
        f.write(response.content)
        data = fits.getdata(f)

    wv = data['wavelength']/10 # convert from Angstroms to nm
    # (erg/cm2/s/Ang)*(1 W/1e7 erg)*(1e3 mW/1 W)*(1e4 cm^2/1 m^2)*(10 Ang/1 nm) = mW/m^2/nm
    F = data['flux_density']*(1/1e7)*(1e3/1)*(1e4/1)*(10/1) # convert from erg/cm2/s/Ang to mW/m^2/nm

    # Remove duplicated wavelengths
    wv, inds = np.unique(wv, return_index=True)
    F = F[inds]

    # Rescale to planet
    F = scale_spectrum_to_planet(wv, F, Teq, stellar_flux)

    # Only consider needed resolution
    if needed_resolution:
        wv, F = rebin_to_needed_resolution(wv, F)

    # Save the spectrum to a file, if desired
    if outputfile is not None:
        save_photochem_spectrum(wv, F, outputfile, scale_to_planet=False)

    return wv, F

def print_hazmat_stars():
    "Prints the stars avaliable in the HAZMAT catalogue"
    print(HAZMAT_STARS_YAML)

HAZMAT_STARS_YAML = \
"""GJ176: {st_logg: 4.78, st_met: 0.147, st_rad: 0.48, st_teff: 3703.49}
GJ436: {st_logg: 4.84, st_met: 0.099, st_rad: 0.42, st_teff: 3586.11}
GJ832: {st_logg: 4.7, st_met: -0.7, st_rad: 0.48, st_teff: 3657.0}
TRAPPIST-1: {st_logg: 5.24, st_met: 0.053, st_rad: 0.12, st_teff: 2566.0}
"""

###
### MUSCLES spectra (https://archive.stsci.edu/prepds/muscles/)
###

def download_muscles_spectrum(star_name, verbose):
    """Downloads a MUSCLES spectrum, then returns the wavelength
    and flux values in photochem units.

    Parameters
    ----------
    star_name : str
        Name of the MUSCLES star.

    Returns
    -------
    wv : ndarray[ndim=1,double]
        Wavelengths of the stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum (mW/m^2/nm)
    """    
    # Look for the spectrum of interest
    all_obs = Observations.query_criteria(
        provenance_name="muscles",
        objectname=star_name,
    )
    data_products = Observations.get_product_list(all_obs)

    # Look through the results for the spectrum we want.
    tmp = None
    for i,a in enumerate(data_products['productFilename']):
        if 'adapt-const-res-sed' in a and '_'+star_name.lower()+'_' in a:
            tmp = data_products[i]
            break

    # If we didn't find it, then we report a problem
    if tmp is None:
        raise Exception('Failed to download '+star_name+' from MUSCLES')

    # Build the URL
    uri = tmp['dataURI']
    base_url = Observations._portal_api_connection.MAST_DOWNLOAD_URL
    url = base_url + "?uri=" + quote(uri, safe=":/")

    if verbose:
        print('Downloading the spectrum at the following URL: '+url)

    # Download
    response = requests.get(url)
    if response.status_code != 200:
        raise Exception('Failed to download '+star_name+' from MUSCLES')

    # Read the download  
    with tempfile.TemporaryFile() as f:
        f.write(response.content)
        data = fits.getdata(f)

    # Get the spectrum
    wv = data['WAVELENGTH']/10 # convert from Angstroms to nm
    # (erg/cm2/s/Ang)*(1 W/1e7 erg)*(1e3 mW/1 W)*(1e4 cm^2/1 m^2)*(10 Ang/1 nm) = mW/m^2/nm
    F = data['FLUX']*(1/1e7)*(1e3/1)*(1e4/1)*(10/1) # convert from erg/cm2/s/Ang to mW/m^2/nm

    # Remove duplicated wavelengths
    wv, inds = np.unique(wv, return_index=True)
    F = F[inds]

    return wv, F

def get_muscles_spectrum(star_name, verbose, nwb=1000):
    """Downloads a MUSCLES spectrum, then adds on a blackbody extending
    the star to 100 microns, and finally rescale the spectrum so that it has
    a total energy consistent with the effective temperature. Returns the 
    wavelength and flux points in photochem units.

    Parameters
    ----------
    star_name : str
        Name of the MUSCLES star
    nwb : int, optional
        Number of blackbody points extending to 100 um, by default 1000

    Returns
    -------
    wv : ndarray[ndim=1,double]
        Wavelengths of the stellar spectrum (nm)
    F : ndarray[ndim=1,double]
        Flux of the stellar spectrum at the stellar surface (mW/m^2/nm)
    """

    # Get the stellar properties
    if star_name not in MUSCLES_STARS:
        raise ValueError(
            'Input `star_name` is not in the MUSCLES library or there is no '+
            'stellar properties for the star. The spelling must match a dictionary'+
            'key in `MUSCLES_STARS`.'
            )
    Teff = MUSCLES_STARS[star_name]['st_teff']

    # Download
    wv, F = download_muscles_spectrum(star_name, verbose)
        
    # Tack on a blackbody to extend the spectrum to 100 um
    wv, F = append_blackbody_to_stellar_spectrum(wv, F, Teff, wv_end=100e3, nwb=nwb)
        
    # Rescale the spectrum so that it's total bolometric flux matches Teff
    factor = stefan_boltzmann(Teff)/energy_in_spectrum(wv, F)
    F *= factor

    return wv, F

def muscles_spectrum(star_name, outputfile=None, Teq=None, stellar_flux=None, needed_resolution=True, verbose=True):
    """Downloads a MUSCLES spectrum (https://archive.stsci.edu/prepds/muscles/), 
    then adds on a blackbody extending the star to 100 microns, and finally 
    rescale the spectrum so that it has a total bolometric insolation at a planet 
    consistent with `Teq` or `stellar_flux`. Finally, the spectrum can be saved to 
    `outputfile` in photochem format.

    Parameters
    ----------
    star_name : str
        Name of the MUSCLES star
    outputfile : str, optional
        If not None, then this is the output filename, be default None.
    Teq : float, optional
        Zero-albedo equilibrium temperature of the planet (K), by default None
    stellar_flux : float, optional
        Stellar flux at the planet (W/m^2), by default None
    needed_resolution : bool, optional
        If True, then the spectrum is rebinned to a resolution 4x higher than
        What is used by the photochemical and climate models which should be
        an adequately high resolution, by default True.
    verbose: bool, optional
        If True, then some information will be printed.
    """
    
    # Download the spectrum
    wv, F = get_muscles_spectrum(star_name, verbose, 1000)

    # Rescale to planet
    F = scale_spectrum_to_planet(wv, F, Teq, stellar_flux)

    # Only consider needed resolution
    if needed_resolution:
        wv, F = rebin_to_needed_resolution(wv, F)

    # Save the spectrum to a file, if desired
    if outputfile is not None:
        save_photochem_spectrum(wv, F, outputfile, scale_to_planet=False)

    return wv, F

def closest_muscles_to_Teff(Teff):
    """Finds the star in the MUSCLES catalogue with an effective temperature
    closest to `Teff`.

    Parameters
    ----------
    Teff : float
        Stellar effective temperature (K)

    Returns
    -------
    dict
        Dictionary with the MUSCLES star name and properties.
    """
    # Get stars with 
    muscles = {a: MUSCLES_STARS[a] for a in MUSCLES_STARS if MUSCLES_STARS[a]['st_teff'] is not None}
    star_names = [key for key in muscles]
    Teffs = np.array([muscles[a]['st_teff'] for a in star_names])
    ind = np.argmin(np.abs(Teff - Teffs))
    star_name = star_names[ind]
    out = muscles[star_name]
    out['name'] = star_name
    return out

def print_muscles_stars():
    "Prints the stars avaliable in the MUSCLES catalogue"
    print(MUSCLES_STARS_YAML)

# All the MUSCLES stars for which we could easily get stellar properties.
# See function `get_muscles_stellar_data` below.
MUSCLES_STARS_YAML = \
"""GJ1132: {st_logg: 5.04, st_met: -0.17, st_rad: 0.22, st_teff: 3229.0}
GJ1214: {st_logg: 5.03, st_met: 0.24, st_rad: 0.22, st_teff: 3101.0}
GJ15A: {st_logg: 4.89, st_met: -0.391, st_rad: 0.38, st_teff: 3742.67}
GJ163: {st_logg: 4.82, st_met: 0.1, st_rad: 0.41, st_teff: 3399.0}
GJ176: {st_logg: 4.78, st_met: 0.147, st_rad: 0.48, st_teff: 3703.49}
GJ436: {st_logg: 4.84, st_met: 0.099, st_rad: 0.42, st_teff: 3586.11}
GJ551: {st_logg: 5.16, st_met: null, st_rad: 0.14, st_teff: 2900.0}
GJ581: {st_logg: 4.94, st_met: -0.088, st_rad: 0.32, st_teff: 3490.39}
GJ649: {st_logg: 4.76, st_met: -0.15, st_rad: 0.5, st_teff: 3734.0}
GJ667C: {st_logg: 4.69, st_met: -0.55, st_rad: null, st_teff: 3650.0}
GJ674: {st_logg: 4.86, st_met: -0.28, st_rad: 0.36, st_teff: 3451.0}
GJ676A: {st_logg: 4.62, st_met: 0.08, st_rad: 0.69, st_teff: 3734.0}
GJ699: {st_logg: 4.9, st_met: -0.56, st_rad: 0.18, st_teff: 3195.0}
GJ729: {st_logg: null, st_met: null, st_rad: null, st_teff: 3240.0}
GJ832: {st_logg: 4.7, st_met: -0.7, st_rad: 0.48, st_teff: 3657.0}
GJ849: {st_logg: 4.8, st_met: 0.09, st_rad: 0.45, st_teff: 3467.0}
GJ876: {st_logg: 4.87, st_met: 0.213, st_rad: 0.37, st_teff: 3293.74}
HAT-P-12: {st_logg: 4.61, st_met: -0.063, st_rad: 0.7, st_teff: 4652.87}
HAT-P-26: {st_logg: 4.5, st_met: 0.028, st_rad: 0.86, st_teff: 5061.9}
HD-149026: {st_logg: 4.17, st_met: 0.33, st_rad: 1.46, st_teff: 6084.0}
HD40307: {st_logg: 4.63, st_met: -0.323, st_rad: 0.72, st_teff: 4867.44}
HD85512: {st_logg: null, st_met: null, st_rad: null, st_teff: 4455.0}
L-678-39: {st_logg: 4.94, st_met: -0.12, st_rad: 0.34, st_teff: 3505.0}
L-98-59: {st_logg: 4.86, st_met: -0.46, st_rad: 0.3, st_teff: 3415.0}
L-980-5: {st_logg: null, st_met: null, st_rad: null, st_teff: 3278.0}
LHS-2686: {st_logg: null, st_met: null, st_rad: null, st_teff: 3220.0}
LP-791-18: {st_logg: 5.12, st_met: -0.09, st_rad: 0.18, st_teff: 2960.0}
TOI-193: {st_logg: 4.47, st_met: 0.25, st_rad: 0.95, st_teff: 5480.0}
TRAPPIST-1: {st_logg: 5.24, st_met: 0.053, st_rad: 0.12, st_teff: 2566.0}
V-EPS-ERI: {st_logg: 4.59, st_met: -0.044, st_rad: 0.76, st_teff: 5020.38}
WASP-127: {st_logg: 4.2, st_met: -0.193, st_rad: 1.35, st_teff: 5828.01}
WASP-17: {st_logg: 4.18, st_met: -0.07, st_rad: 1.57, st_teff: 6548.0}
WASP-43: {st_logg: 4.5, st_met: -0.01, st_rad: 0.65, st_teff: 4500.0}
WASP-77A: {st_logg: 4.48, st_met: -0.1, st_rad: 0.91, st_teff: 5617.0}
"""
MUSCLES_STARS = yaml.safe_load(MUSCLES_STARS_YAML)

# The function below is only needed when updating this file.

def get_muscles_stellar_data():
    """Queries the Nasa Exoplant Archive to try to get stellar properties for 
    each the MUSCLES stars. The function outputs a yaml file for all stars for which 
    properties were found.
    """
    from bs4 import BeautifulSoup
    from ._format import flowmap

    # All the MUSCLES stars, with some alternative names
    ALL_MUSCLES_STARS = {
        'GJ1132': [],
        'GJ1214': [],
        'GJ15A': [],
        'GJ163': [],
        'GJ176': ['Gaia DR2 3409711211681795584'],
        'GJ436': [],
        'GJ581': [],
        'GJ649': [],
        'GJ667C': [],
        'GJ674': [],
        'GJ676A': [],
        'GJ699': ['Gaia DR2 4472832130942575872'],
        'GJ729': ['Gaia DR2 4075141768785646848',],
        'GJ832': [],
        'GJ849': [],
        'GJ876': [],
        'HD40307': [],
        'HD85512': ['Gaia DR2 5412947081287925504'],
        'L-980-5': ['Gaia DR2 3652796572020424448'],
        'LHS-2686': [],
        'TRAPPIST-1': [],
        'V-EPS-ERI': ['Gaia DR2 5164707970261630080'],
        'GJ551': ['Gaia DR2 5853498713160606720'],
        'LP-791-18': ['LP791-18'],
        'WASP-43': [],
        'HD-149026': ['HD149026'],
        'L-98-59': ['L98-59'],
        'HAT-P-26': [],
        'L-678-39': ['Gaia DR2 5664814198431308288'],
        'TOI-193': ['Gaia DR2 2307504062045842688'],
        'WASP-17': [],
        'WASP-127': [],
        'HAT-P-12': [],
        'WASP-77A': [],
    }

    # NEA api
    response = requests.get(
        "https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query="+
        "select+sy_name,hostname,gaia_id,st_teff,st_met,st_rad,st_logg,st_refname+from+stellarhosts&format=json"
    )
    data = response.json()

    # First find stars for which there is data in the NEA
    exo_names = np.array([a['hostname'].replace(' ','') for a in data])
    gaia_names = np.array([a['gaia_id'] for a in data])
    stars = {}
    for key in ALL_MUSCLES_STARS:
        # Search using primary name
        inds = np.where(key == exo_names)[0]
        if len(inds) > 0:
            stars[key] = inds
            continue
        
        # Search using alternative name
        inds = np.array([],dtype=np.int64)
        for altname in ALL_MUSCLES_STARS[key]:
            inds = np.where(altname == exo_names)[0]
            if len(inds) > 0:
                break
            inds = np.where(altname == gaia_names)[0]
            if len(inds) > 0:
                break
        if len(inds) > 0:
            stars[key] = inds
            continue

        # The star was not found
        stars[key] = []

    # Next, sort all the data by publication year (newest to oldest)
    for key in stars:
        years = []
        for ind in stars[key]:
            soup = BeautifulSoup(data[ind]['st_refname'],features="html.parser")
            ref = soup.find(href=True)['refstr']
            year = 0
            try:
                year = int(ref.split('__')[-1])
            except ValueError:
                pass
            years.append(year)
        if len(stars[key]) > 0:
            inds1 = np.argsort(years)[::-1]
            stars[key] = list(np.array(stars[key])[inds1])

    # Next, sort so that publications will all data of interest are first
    datakeys = ['st_teff','st_met','st_rad','st_logg']
    for key in stars:
        best_inds = []
        other_inds = []
        for ind in stars[key]:
            all_info = True
            for key1 in datakeys:
                if data[ind][key1] is None:
                    all_info = False
                    break
            if all_info:
                best_inds.append(ind)
            else:
                other_inds.append(ind)
        if len(stars[key]) > 0:
            stars[key] = best_inds + other_inds

    # Extract the stellar properties
    stars_new = {}
    for key in stars:
        res = {key1: None for key1 in datakeys}
        for ind in stars[key]:
            for key1 in datakeys:
                if res[key1] is None and data[ind][key1] is not None:
                    res[key1] = data[ind][key1]
        stars_new[key] = res
    stars = stars_new

    # Save as a yaml file
    for star in stars:
        stars[star] = flowmap(stars[star])
    with open('muscles_stars.yaml','w') as f:
        yaml.dump(stars,f,yaml.Dumper)
