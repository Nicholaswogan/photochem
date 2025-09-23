import numpy as np
import numba as nb
from numba import types
import yaml
from tempfile import NamedTemporaryFile
import sys

from scipy import constants as const
from astropy import constants
import astropy.units as u
from stsynphot import grid_to_spec
from synphot import units

from .._clima import AdiabatClimate, ClimaException, rebin, rebin_with_errors
from .. import utils
from ..utils import stars

# Optional imports
try:
    from picaso import justdoit as jdi
    import pandas as pd
except ModuleNotFoundError:
    pass
try:
    from pandexo.engine import justdoit as pandexo_jdi
except ModuleNotFoundError:
    pass

# We use a jitclass because it has strict type checking
@nb.experimental.jitclass()
class ThermalEmissionData():
    Teq : types.double # Zero albedo equilibrium temperature # type: ignore
    R_planet : types.double # Radius of planet in cm # type: ignore
    M_planet : types.double # Mass of planet in grams # type: ignore
    R_star : types.double # Radius of star in cm # type: ignore
    wavl_star : types.double[:] # Wavelength grid for stellar surface (um) # type: ignore
    flux_star : types.double[:] # Flux of stellar surface (ergs/cm^2/s/cm) # type: ignore
    flux_star_c : types.double[:] # Flux of stellar surface regrided to Clima IR bins (ergs/cm^2/s/cm) # type: ignore
    def __init__(self):
        pass

class AdiabatClimateThermalEmission(AdiabatClimate):
    "An extension of the AdiabatClimate class to interpret thermal emission observations."

    def __init__(self, Teq, M_planet, R_planet, R_star, Teff=None, metal=None, logg=None, 
                 stellar_surface_file=None, catdir='phoenix', stellar_surface_scaling=1.0,
                 species=None, condensates=None, species_file=None, 
                 opacities_file=None, data_dir=None, nz=50, number_of_zeniths=4):
        """Initializes the code. 

        Parameters
        ----------
        Teq : float
            Zero albedo equilibrium temperature of the planet (K).
        M_planet : float
            Mass of the planet in Earth masses.
        R_planet : float
            Radius of the planet in Earth radii.
        R_star : float
            Stellar radius in solar radii.
        Teff : float, optional
            Stellar effective temperature in K.
        metal : float, optional
            log10 metallicity of the star.
        logg : float, optional
            log10 gravity of the star in cgs units.
        stellar_surface_file : str, optional
            Path to the stellar surface file. If this file is not specified, then the code will
            alternatively construct a stellar spectrum using pysynphot. The file must have two 
            columns where the first is wavelength in nm, and the second is stellar flux in mW/m^2/nm. 
            The first line of the file is always skipped, with the assumption they are column labels.
        catdir : str, optional
            The stellar database, by default 'phoenix'
        stellar_surface_scaling : float, optional
            Optional scaling to apply to the stellar surface flux, by default 1.0
        species_file : str, optional
            Path to a settings file. If None, then a default file is used.
        opacities_file : str, optional
            Path to a file describing the used opacities. If None, then a default file is used.
        data_dir : str, optional
            Path to where climate model data is stored. If None, then installed data is used.
        nz : int, optional
            Number of vertical layers in the climate model, by default 50
        number_of_zeniths : int, optional
            Number of zenith angles in the radiative transfer calculation, by default 1
        """ 
        
        # Settings file
        if species is None:
            species = ['H2O', 'CO2', 'O2', 'SO2']
        if condensates is None:
            condensates = ['H2O', 'CO2']
        if species_file is None:
            species_dict = utils.species_dict_for_climate(species, condensates)
        else:
            with open(species_file,'r') as f:
                species_dict = yaml.load(f, Loader=yaml.Loader)

        # Settings file
        opacities = {
            'k-distributions': True, 
            'CIA': True, 
            'rayleigh': True, 
            'photolysis-xs': True
        }
        settings_dict = utils.settings_dict_for_climate(
            planet_mass=float(M_planet*constants.M_earth.to('g').value), 
            planet_radius=float(R_planet*constants.R_earth.to('cm').value), 
            surface_albedo=0.0, 
            number_of_layers=int(nz), 
            number_of_zenith_angles=int(number_of_zeniths), 
            photon_scale_factor=1.0, 
            opacities=opacities
        )
        if opacities_file is not None:
            # If custom optical properties are specified, then we use them
            with open(opacities_file,'r') as f:
                optical_properties = yaml.load(f, Loader=yaml.Loader)
            settings_dict['optical-properties'] = optical_properties['optical-properties']

        # Stellar flux
        if stellar_surface_file is None:
            if Teff is None or metal is None or logg is None:
                raise ClimaException('If `stellar_surface_file` is None, then `Teff`, `metal` and `logg` '+
                                     'must all be supplied')
            # Get stellar flux from pysynphot
            wv_star, F_star, wv_planet, F_planet = make_pysynphot_stellar_spectrum(Teq, Teff, metal, logg, catdir)
        else:
            # Get stellar flux from a file
            wv_star, F_star, wv_planet, F_planet = make_file_stellar_spectrum(Teq, stellar_surface_file)

        F_star = F_star*stellar_surface_scaling # Scale stellar surface, if necessary
        # Load the flux at the planet to a string
        flux_str = stars.photochem_spectrum_string(wv_planet, F_planet, scale_to_planet=False)

        # Initialize
        with NamedTemporaryFile('w') as f_species:
            # Write species file
            yaml.safe_dump(species_dict, f_species)
            with NamedTemporaryFile('w') as f_settings:
                # Write settings file
                yaml.safe_dump(settings_dict, f_settings)
                with NamedTemporaryFile('w') as f_flux:
                    # Write stellar flux file
                    f_flux.write(flux_str)
                    f_flux.flush()
                    
                    # Initialize AdiabatClimate
                    super().__init__(
                        f_species.name, 
                        f_settings.name, 
                        f_flux.name,
                        data_dir=data_dir
                    )

        # Change default parameters
        self.solve_for_T_trop = True # Enable solving for T_trop
        self.tidally_locked_dayside = True # Enable Tidally locked dayside calculations
        self.max_rc_iters = 30 # Lots of iterations
        self.P_top = 10.0 # 10 dynes/cm^2 top, or 1e-5 bars.

        # Tack on new object
        self.thermdat = ThermalEmissionData()

        # Save some information
        self.thermdat.Teq = Teq
        self.thermdat.R_planet = float(R_planet*constants.R_earth.to('cm').value) # Planet radius in cm
        self.thermdat.M_planet = float(M_planet*constants.M_earth.to('g').value) # Planet mass in grams
        self.thermdat.R_star = float(R_star*constants.R_sun.to('cm').value) # Star Radius in cm

        # Rebin stellar flux to Clima IR grid
        self.thermdat.wavl_star = stars.make_bins(wv_star) # Stellar wavelength grid (microns)
        self.thermdat.flux_star = F_star # Stellar flux in each bin (ergs/cm^2/s/cm)
        self.thermdat.flux_star_c = rebin(self.thermdat.wavl_star, self.thermdat.flux_star, self.rad.ir.wavl/1e3)

        # Placeholder for picaso object
        self.ptherm = None
    
    def surface_temperature_robust(self, P_i, T_guess_mid=None, T_perturbs=None):
        """Similar to self.surface_temperature, except more robust.

        Parameters
        ----------
        P_i : ndarray[double,ndim=1]
            Array of surface pressures of each species (dynes/cm^2)
        T_guess_mid : float, optional
            Middle T_surface guess, by default None
        T_perturbs : float, optional
            Perturbations to `T_guess_mid` to attempt, by default None

        Returns
        -------
        bool
            If True, the the model is converged.
        """        

        if T_guess_mid is None:
            T_guess_mid = self.thermdat.Teq*1.5
        
        if T_perturbs is None:
            T_perturbs = np.array([0.0, 50.0, -50.0, 100.0, -100.0, 150.0, 800.0, 600.0, 400.0, 300.0, 200.0])

        for i,T_perturb in enumerate(T_perturbs):
            T_surf_guess = T_guess_mid + T_perturb
            try:
                self.T_trop = self.rad.skin_temperature(0.0)*1.2
                self.surface_temperature(P_i, T_surf_guess)
                converged = True
                break
            except ClimaException as e:
                converged = False
        
        return converged

    def RCE_robust(self, P_i, T_guess_mid=None, T_perturbs=None):
        """Similar to self.RCE, except more robust.

        Parameters
        ----------
        P_i : ndarray[double,ndim=1]
            Array of surface pressures of each species (dynes/cm^2)
        T_guess_mid : float, optional
            Middle temperature guess, by default None
        T_perturbs : float, optional
            Perturbations to `T_guess_mid` to attempt, by default None

        Returns
        -------
        bool
            If True, the the model is converged.
        """    

        if T_guess_mid is None:
            T_guess_mid = self.thermdat.Teq*1.5
        
        if T_perturbs is None:
            T_perturbs = np.array([0.0, 50.0, -50.0, 100.0, -100.0, 150.0, 800.0, 600.0, 400.0, 300.0, 200.0])

        # First, we try a bunch of isothermal atmospheres.
        for i,T_perturb in enumerate(T_perturbs):
            T_surf_guess = T_guess_mid + T_perturb
            T_guess = np.ones(self.T.shape[0])*T_surf_guess
            try:
                converged = self.RCE(P_i, T_surf_guess, T_guess)
                if converged:
                    break
            except ClimaException:
                converged = False

        # If we success, then we return.
        if converged:
            return converged

        # If not converged, then we run the simple climate model.
        converged_simple = self.surface_temperature_robust(P_i)
        if not converged_simple:
            # If this fails, then we give up, returning no convergence
            return False

        # If simple climate model converged, then save the atmosphere
        T_surf_guess, T_guess, convecting_with_below_guess = self.T_surf, self.T, self.convecting_with_below

        # Initial guess without convecting pattern
        try:
            converged = self.RCE(P_i, T_surf_guess, T_guess)
        except ClimaException:
            converged = False

        # If we success, then we return.
        if converged:
            return converged
        
        # If not converged, then we try again with the convecting pattern
        try:
            converged = self.RCE(P_i, T_surf_guess, T_guess, convecting_with_below_guess)
        except ClimaException:
            converged = False

         # Success or not, we return at this point.
        return converged

    def set_custom_albedo(self, wv, albedo):
        """Sets a cutsom surface albedo/emissivity. The input is 
        constantly extrapolated.

        Parameters
        ----------
        wv : ndarray[double,ndim=1]
            Wavelength points in microns.
        albedo : ndarray[double,ndim=1]
            Surface albedo at each wavelength.
        """        

        # Get wv grid of IR
        freq = self.rad.ir.freq # Hz
        freq_av = (freq[1:]+freq[:-1])/2 # Hz
        wv_ir_av = 1e6*const.c/freq_av # microns

        # Get wv grid of Solar
        freq = self.rad.sol.freq # Hz
        freq_av = (freq[1:]+freq[:-1])/2 # Hz
        wv_sol_av = 1e6*const.c/freq_av # microns

        albedo_c = np.interp(wv_sol_av, wv, albedo)
        emissivity_c = 1.0 - np.interp(wv_ir_av, wv, albedo)

        self.rad.surface_albedo = albedo_c
        self.rad.surface_emissivity = emissivity_c

        # Also set PICASO surface albedo
        if self.ptherm is not None:
            self.set_custom_albedo_picaso(wv, albedo)

###
### Making spectra
###

    def fpfs_blackbody(self, wavl, T, albedo):
        """Generates the planet-to-star flux ratio for a blackbody of a given
        temperature.

        Parameters
        ----------
        wavl : ndarray[double,ndim=1]
            Edges of the wavelength grid in microns.
        T : float
            Blackbody temperature in K.
        albedo : float
            Albedo of the surface.

        Returns
        -------
        F_planet : ndarray[double,ndim=1]
            Flux of the planet in each wavelength bin (ergs/cm^2/s/cm)
        F_star : ndarray[double,ndim=1]
            Flux of the star in each wavelength bin (ergs/cm^2/s/cm)
        fpfs : ndarray[double,ndim=1]
            Planet-to-star flux ratio in each wavelength bin.
        """        
        thermdat = self.thermdat

        wv_av = (wavl[1:] + wavl[:-1])/2
        F_planet = (1 - albedo)*stars.blackbody_cgs(T, wv_av/1e4)*np.pi
        F_star = rebin(thermdat.wavl_star, thermdat.flux_star, wavl)
        fpfs = F_planet/F_star * (thermdat.R_planet**2/thermdat.R_star**2)

        return F_planet, F_star, fpfs
    
    def fpfs_instant_reradiation(self, wavl, albedo):
        flux = stars.equilibrium_temperature_inverse(self.thermdat.Teq, albedo)
        T = bare_rock_dayside_temperature(flux, albedo, 2/3)
        return self.fpfs_blackbody(wavl, T, albedo)

    def fpfs(self):
        """Generates the planet-to-star flux ratio using the most recent
        radiative transfer results.

        Returns
        -------
        wavl : ndarray[double,ndim=1]
            Edges of the wavelength grid in microns.
        F_planet : ndarray[double,ndim=1]
            Flux of the planet in each wavelength bin (ergs/cm^2/s/cm)
        fpfs : ndarray[double,ndim=1]
            Planet-to-star flux ratio in each wavelength bin.
        """  

        freq = self.rad.ir.freq # Hz
        freq_av = (freq[1:]+freq[:-1])/2 # Hz
        wv_av = 1e6*const.c/freq_av # microns
        # c.rad.wrk_ir.fup_a is mW/m^2/Hz. Here I convert to W/m^2/um
        F1 = 1e-3*self.rad.wrk_ir.fup_a[-1,:]*(freq_av/(wv_av))
        F_planet = F1*(1e7/1)*(1/1e4)*(1e6/1)*(1/1e2) # convert to erg/cm^2/s/cm

        fpfs_c = F_planet/self.thermdat.flux_star_c * (self.thermdat.R_planet**2/self.thermdat.R_star**2)

        return self.rad.ir.wavl/1e3, F_planet, fpfs_c
    
    def brightness_temperature(self):
        """Returns the brightness temperature using the most recent radiative transfer
        result

        Returns
        -------
        wavl : ndarray[double,ndim=1]
            Edges of the wavelength grid in microns.
        T : ndarray[double,ndim=1]
            The brightness temperature (K) in each wavelength bin.
        """        
        wavl, F_planet, _ = self.fpfs()
        T = np.empty(F_planet.shape[0])
        for i in range(F_planet.shape[0]):
            wv_av = (wavl[i] + wavl[i+1])/2
            T[i] = inverse_blackbody(F_planet[i]/np.pi,np.array([wv_av/1e4]))
        return wavl, T

###
### PICASO utilities
###

    def initialize_picaso_from_clima(self, filename_db, opannection_kwargs={}):

        if 'picaso' not in sys.modules:
            raise Exception('To use picaso, you must first install it.')

        # Planet and star
        M_planet = self.thermdat.M_planet/constants.M_earth.to('g').value
        R_planet = self.thermdat.R_planet/constants.R_earth.to('cm').value
        R_star = self.thermdat.R_star/constants.R_sun.to('cm').value

        # Stellar flux
        # ergs/cm^2/s/cm. Must convert to erg/cm^2/s/A (FLAM units)
        # ergs/cm^2/s/cm * (1e2 cm / m) * (m / 1e10 A)
        flux_star = self.thermdat.flux_star*(1e2/1)*(1/1e10)
        wv_star = (self.thermdat.wavl_star[1:] + self.thermdat.wavl_star[:-1])/2 # microns

        with NamedTemporaryFile('w') as f:
            np.savetxt(f, np.concatenate(([wv_star], [flux_star])).T)
            f.flush()

            star_kwargs = {
                'filename': f.name,
                'w_unit': 'um',
                'f_unit': 'FLAM'
            }
            self.ptherm = PicasoThermalEmission(
                filename_db, 
                M_planet, 
                R_planet, 
                R_star, 
                opannection_kwargs=opannection_kwargs,
                star_kwargs=star_kwargs
            )

        return self.ptherm
    
    def make_picaso_atm(self):
        atm = {}
        atm['pressure'] = np.append(self.P_surf,self.P)/1e6
        atm['temperature'] = np.append(self.T_surf,self.T)
        for i,sp in enumerate(self.species_names):
            atm[sp] = self.f_i[:,i]
            atm[sp] = np.append(self.f_i[0,i],atm[sp])
        for key in atm:
            atm[key] = atm[key][::-1].copy()
        atm = pd.DataFrame(atm)
        return atm
    
    def set_custom_albedo_picaso(self, wv, albedo):
        if self.ptherm is None:
            raise Exception('You must first initialize picaso with `initialize_picaso_from_clima`')
        self.ptherm.set_custom_albedo(wv, albedo)

    def fpfs_picaso(self, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        if self.ptherm is None:
            raise Exception('You must first initialize picaso with `initialize_picaso_from_clima`')
        return self.ptherm.fpfs(self.make_picaso_atm(), R=R, wavl=wavl, atmosphere_kwargs=atmosphere_kwargs, **kwargs)
    
    def brightness_temperature_picaso(self, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        if self.ptherm is None:
            raise Exception('You must first initialize picaso with `initialize_picaso_from_clima`')
        return self.ptherm.brightness_temperature(self.make_picaso_atm(), R=R, wavl=wavl, atmosphere_kwargs=atmosphere_kwargs, **kwargs)

    def rprs2_picaso(self, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        if self.ptherm is None:
            raise Exception('You must first initialize picaso with `initialize_picaso_from_clima`')
        return self.ptherm.rprs2(self.make_picaso_atm(), R=R, wavl=wavl, atmosphere_kwargs=atmosphere_kwargs, **kwargs)

###
### PandExo utilities
###

    def create_exo_dict(self, total_observing_time, eclipse_duration, kmag, starpath, calculation):

        if 'pandexo' not in sys.modules:
            raise Exception('To use pandexo, you must first install it.')

        exo_dict = pandexo_jdi.load_exo_dict()

        exo_dict['observation']['sat_level'] = 80
        exo_dict['observation']['sat_unit'] = '%'
        exo_dict['observation']['noccultations'] = 1
        exo_dict['observation']['R'] = None
        exo_dict['observation']['baseline_unit'] = 'total'
        exo_dict['observation']['baseline'] = total_observing_time
        exo_dict['observation']['noise_floor'] = 0

        exo_dict['star']['type'] = 'user'
        exo_dict['star']['mag'] = kmag
        exo_dict['star']['ref_wave'] = 2.22
        exo_dict['star']['starpath'] = starpath
        exo_dict['star']['w_unit'] = 'um'
        exo_dict['star']['f_unit'] = 'erg/cm2/s/Hz'
        exo_dict['star']['radius'] = self.thermdat.R_star
        exo_dict['star']['r_unit'] = 'cm'

        exo_dict['planet']['type'] = 'constant'
        exo_dict['planet']['transit_duration'] = eclipse_duration
        exo_dict['planet']['td_unit'] = 's'
        exo_dict['planet']['radius'] = self.thermdat.R_planet
        exo_dict['planet']['r_unit'] = 'cm'
        if calculation == 'thermal':
            exo_dict['planet']['f_unit'] = 'fp/f*'
            stellar_flux = stars.equilibrium_temperature_inverse(self.thermdat.Teq, 0.0)
            T_day = bare_rock_dayside_temperature(stellar_flux, 0.0, 2/3)
            exo_dict['planet']['temp'] = T_day
        elif calculation == 'transmission':
            exo_dict['planet']['f_unit'] = 'rp^2/r*^2'
        else:
            raise Exception('`calculation` must be "thermal" or "transmission"')

        return exo_dict
    
    def pandexo_stellar_filestr(self):

        wavl = self.thermdat.wavl_star # microns
        wv_av = (wavl[1:] + wavl[:-1])/2 
        freq_av = const.c/(wv_av*1e-6)
        F = self.thermdat.flux_star # ergs/cm^2/s/cm
        # ergs/cm^2/s/cm * (1e2 cm / 1 m) * (1 m / 1e6 um) = ergs/cm^2/s/um
        F = F*(1e2/1)*(1/1e6)
        # ergs/cm^2/s/um * (um/Hz) = ergs/cm^2/s/Hz
        F = F*(wv_av/freq_av)

        filestr = ''
        fmt = '{:20}'
        for i in range(wv_av.shape[0]):
            filestr += fmt.format('%e'%wv_av[i])+fmt.format('%e'%F[i])+'\n'

        return filestr
    
    def _run_pandexo(self, total_observing_time, eclipse_duration, kmag, inst, calculation, verbose=False, **kwargs):

        if 'pandexo' not in sys.modules:
            raise Exception('To use pandexo, you must first install it.')

        with NamedTemporaryFile('w') as f:
            
            # Write the Stellar flux to a file
            f.write(self.pandexo_stellar_filestr())
            f.flush()

            # Create the exo dict
            exo_dict = self.create_exo_dict(total_observing_time, eclipse_duration, kmag, f.name, calculation)

            # Run pandexo
            result = pandexo_jdi.run_pandexo(exo_dict, inst, verbose=verbose, **kwargs)

        return result

    def run_pandexo(self, total_observing_time, eclipse_duration, kmag, inst, calculation='thermal', R=None, ntrans=1, verbose=False, **kwargs):

        # inst is just a string
        assert isinstance(inst, str)
        result = self._run_pandexo(total_observing_time, eclipse_duration, kmag, [inst], calculation, verbose, **kwargs)

        spec = result['FinalSpectrum']
        wavl = stars.make_bins(spec['wave'])
        F = spec['spectrum']
        err = spec['error_w_floor']
        err = err/np.sqrt(ntrans)

        if R is not None:
            wavl_n = stars.grid_at_resolution(np.min(wavl), np.max(wavl), R)
            F_n, err_n = rebin_with_errors(wavl, F, err, wavl_n)
            wavl = wavl_n
            F = F_n
            err = err_n

        return wavl, F, err, result

###
### Making stellar fluxes
###

def get_planet_flux(Teq, wv_star, F_star):

    # Convert to units in climate model
    wv_0 = wv_star*1e3 # to nm
    # erg/cm^2/s/cm * (W/(erg/s)) * (mW/W) * (cm^2/m^2) * (cm/m) * (m/nm) = mW/m^2/nm
    F_0 = F_star*(1/1e7)*(1e3/1)*(1e4/1)*(1e2/1)*(1/1e9) 

    # Interpolate to smaller resolution appropriate for climate modeling
    wv_planet, F_planet = stars.rebin_to_needed_resolution(wv_0, F_0)

    # Rescale
    F_planet = stars.scale_spectrum_to_planet(wv_planet, F_planet, Teq=Teq)
    
    return wv_planet, F_planet

def make_pysynphot_stellar_spectrum(Teq, Teff, metal, logg, catdir='phoenix'):
    """Create stellar spectrum for AdiabatClimate and thermal emission predictions
    using the pysynphot package

    Parameters
    ----------
    Teq : float
        Planet equilibrium temperature assuming zero bond albedo.
    Teff : float
        Stellar effective temperature
    metal : float
        Stellar metallicity
    logg : float
        Stellar gravity in log space.
    catdir : str, optional
        Stellar database, by default 'phoenix'

    Returns
    -------
    tuple
        `(wv_star, F_star, wv_planet, F_planet)`, where wv_star is the wavelengths
        of the stellar flux in microns, F_star is the stellar surface flux in ergs/cm^2/s/cm,
        wv_planet is the wavelengths of the flux at the planet in nm, and 
        F_planet is the stellar flux at the planet in mW/m^2/nm.
    """    
    # Get spectrum
    sp = grid_to_spec(catdir, Teff, metal, logg)
    wv_star = sp.waveset.to(u.um).value # um
    F_star = units.convert_flux(sp.waveset, sp(sp.waveset), units.FLAM).value*1e8 # Convert to ergs/cm2/s/cm
    
    wv_planet, F_planet = get_planet_flux(Teq, wv_star, F_star)

    return wv_star, F_star, wv_planet, F_planet

def make_file_stellar_spectrum(Teq, stellar_surface_file):

    # Get the spectrum. wv_star is nm and F_star is mW/m^2/nm. 
    # Must convert to microns and ergs/cm2/s/cm
    wv_star, F_star = np.loadtxt(stellar_surface_file, skiprows=1).T

    # nm * (1 um/ 1e3 nm) = um
    wv_star = wv_star/1e3

    # mW/m^2/nm * (1e9 nm / m) * (m / 1e2 cm) * (W / 1e3 mW) * (1e7 (ergs/s) / W) * (m^2 / 1e4 cm^2)
    # = ergs/cm^2/s/cm
    F_star = F_star*(1e9/1)*(1/1e2)*(1/1e3)*(1e7/1)*(1/1e4)

    wv_planet, F_planet = get_planet_flux(Teq, wv_star, F_star)

    return wv_star, F_star, wv_planet, F_planet

@nb.njit(types.double[:](types.double,types.double[:]))
def inverse_blackbody(B, lam):
    """
    Inverse of the `blackbody` function.

    Parameters
    ----------
    B : float
        Blackbody flux in erg/cm^2/s/cm/sr
    lam : ndarray[ndim=1,double]
        Wavelength (cm)
    
    Returns
    -------
    ndarray[ndim=1,double]
        The temperature at the input wavelengths in K.
    """
    h = const.h*1e7 # erg s 
    c = const.c*1e2 # cm/s
    k = const.k*1e7 # erg / K
    T = ((c*h)/(k*lam))*(1/np.log(1 + (2*c**2*h)/(B*lam**5)))
    return T

@nb.njit()
def bare_rock_dayside_temperature(stellar_radiation, bond_albedo, f_term):
    T_eq = stars.equilibrium_temperature(stellar_radiation, bond_albedo)
    return T_eq*(4*f_term)**(1/4) 

class PicasoThermalEmission():

    def __init__(self, filename_db, M_planet, R_planet, R_star, opannection_kwargs={}, star_kwargs={}):

        self.opa = jdi.opannection(filename_db=filename_db, **opannection_kwargs)
        self.case = jdi.inputs()
        self.case.phase_angle(0)
        self.case.gravity(mass=M_planet, mass_unit=jdi.u.Unit('M_earth'),
                     radius=R_planet, radius_unit=jdi.u.Unit('R_earth'))
        self.case.star(self.opa, radius=R_star, radius_unit=jdi.u.Unit('R_sun'), **star_kwargs)
        self.case.surface_reflect(np.ones(self.opa.wno.shape[0])*0.0,self.opa.wno)

    def set_custom_albedo(self, wv, albedo):
        self.case.surface_reflect(albedo[::-1].copy(), self.opa.wno, (1e4/wv[::-1]).copy())

    def _spectrum(self, atm, calculation='thermal', atmosphere_kwargs={}, **kwargs):
        self.case.atmosphere(df=atm, verbose=False, **atmosphere_kwargs)
        self.case.approx(p_reference=np.max(atm['pressure'].to_numpy()))
        df = self.case.spectrum(self.opa, calculation=calculation, **kwargs)
        return df

    def fpfs(self, atm, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        df = self._spectrum(atm, atmosphere_kwargs=atmosphere_kwargs, **kwargs)

        wavl_h = stars.make_bins(1e4/df['wavenumber'][::-1].copy())
        fpfs_h = df['fpfs_thermal'][::-1].copy()
        fp_h = df['thermal'][::-1].copy()

        if wavl is None:
            wavl = stars.grid_at_resolution(np.min(wavl_h), np.max(wavl_h), R)

        fp = rebin(wavl_h, fp_h, wavl)
        fpfs = rebin(wavl_h, fpfs_h, wavl)

        return wavl, fp, fpfs
    
    def brightness_temperature(self, atm, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        wavl, F_planet, _ = self.fpfs(atm, R=R, wavl=wavl, atmosphere_kwargs=atmosphere_kwargs, **kwargs)
        T = np.empty(F_planet.shape[0])
        for i in range(F_planet.shape[0]):
            wv_av = (wavl[i] + wavl[i+1])/2
            T[i] = inverse_blackbody(F_planet[i]/np.pi,np.array([wv_av/1e4]))
        return wavl, T
    
    def rprs2(self, atm, R=100, wavl=None, atmosphere_kwargs={}, **kwargs):
        df = self._spectrum(atm, 'transmission', atmosphere_kwargs=atmosphere_kwargs, **kwargs)

        wavl_h = stars.make_bins(1e4/df['wavenumber'][::-1].copy())
        rprs2_h = df['transit_depth'][::-1].copy()

        if wavl is None:
            wavl = stars.grid_at_resolution(np.min(wavl_h), np.max(wavl_h), R)

        rprs2 = rebin(wavl_h, rprs2_h, wavl)

        return wavl, rprs2
