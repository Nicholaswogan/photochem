import numpy as np

from tempfile import NamedTemporaryFile

from photochem.utils import climate
from photochem.utils._format import yaml
from photochem.clima import AdiabatClimate, ClimaException
from photochem import EvoAtmosphere, PhotoException

class CoupledPhotochemistryClimate():

    def __init__(
            self, 
            mechanism_file, 
            settings_file, 
            flux_file, 
            climate_condensates=('H2O', 'CO2'), 
            climate_number_of_layers=50,
            climate_top_pressure=1.0,
            data_dir=None,
        ):

        # Load files for later use
        with open(mechanism_file, 'r') as f:
            self.mechanism = yaml.load(f, yaml.Loader)
        with open(settings_file, 'r') as f:
            self.settings = yaml.load(f, yaml.Loader)

        # Initialize photochemistry code
        self.pc = EvoAtmosphere(
            mechanism_file=mechanism_file, 
            settings_file=settings_file, 
            flux_file=flux_file, 
            data_dir=data_dir
        )

        # Extract information from photochemistry species file
        climate_species = []
        for i,sp in enumerate(self.mechanism['species']):
            climate_species.append(sp['name'])

        # Validate and save the climate condensibles
        climate_condensates = tuple(climate_condensates)
        self._save_climate_condensates(climate_condensates)

        with NamedTemporaryFile(suffix='.yaml') as f1:
            climate.species_file_for_climate(
                input_filename=mechanism_file,
                filename=f1.name, 
                species=climate_species, 
                condensates=climate_condensates
            )
            f1.flush()
            with NamedTemporaryFile(suffix='.yaml') as f2:
                climate.settings_file_for_climate(
                    filename=f2.name, 
                    planet_mass=float(self.settings['planet']['planet-mass']), 
                    planet_radius=float(self.settings['planet']['planet-radius']), 
                    surface_albedo=float(self.settings['planet']['surface-albedo']), 
                    number_of_layers=climate_number_of_layers,
                    number_of_zenith_angles=4, 
                    photon_scale_factor=float(self.settings['planet'].get('photon-scale-factor', 1.0))
                )
                f2.flush()

                self.c = AdiabatClimateForCoupled(
                    species_file=f1.name, 
                    settings_file=f2.name, 
                    flux_file=flux_file,
                    top_pressure=climate_top_pressure, 
                    data_dir=data_dir,
                )

        # Make the RH the same between models
        self._set_climate_rh_from_photochem()

    def _set_climate_rh_from_photochem(self):

        particles = self.mechanism.get('particles', [])
        gas_phase = {}
        for i, particle in enumerate(particles):
            if particle['formation'] == 'saturation':
                gas_phase[particle['name']] = particle['gas-phase']
        RH = np.ones_like(self.c.RH)
        for i,sp in enumerate(self.pc.dat.species_names[:self.pc.dat.np]):
            sp1 = gas_phase.get(sp)
            if sp1 in self.c.species_names:
                ind = self.c.species_names.index(sp1)
                RH[ind] = self.pc.var.cond_params[i].RHc
        self.c.RH = RH

    def _save_climate_condensates(self, climate_condensates):

        photochem_condensates = {
            particle['gas-phase']
            for particle in self.mechanism.get('particles', [])
            if particle['formation'] == 'saturation'
        }
        unsupported = set(climate_condensates) - photochem_condensates
        if unsupported:
            raise ValueError(
                'Each climate condensate must have a saturation particle in '
                'the photochemical mechanism. Missing for: '
                + ', '.join(sorted(unsupported))
            )

        self.climate_condensates = climate_condensates

    def solution(self):
        """Return the current photochemical atmosphere without modifying it.

        This does not check chemical or coupled convergence.

        Returns
        -------
        dict[str, numpy.ndarray]
            Layer-center altitude (``alt``, cm), temperature (``temp``, K),
            pressure (dyn/cm²), gas density (molecules/cm³), and species
            abundances, following ``pc.mole_fraction_dict`` conventions.
        """
        if not self.pc.atmosphere_initialized:
            raise PhotoException('The photochemical atmosphere has not been initialized.')
        return self.pc.mole_fraction_dict()

    def profiles_from_climate(self):
        """Return bottom-to-top pressure (dyn/cm²), temperature (K), and gas VMRs.

        Profiles include an explicit surface point followed by climate layers.
        Returned arrays are independent of the climate object's state.
        """
        c = self.c
        pressure = np.append(c.P_surf, c.P)
        temperature = np.append(c.T_surf, c.T)
        f_i = np.concatenate((c.f_i_surf.reshape(1,len(c.f_i_surf)),c.f_i))
        mix = {}
        for i,sp in enumerate(c.species_names):
            mix[sp] = f_i[:,i]
        return pressure, temperature, mix

    def profiles_from_photochem(self):
        """Return bottom-to-top layer pressure (dyn/cm²), temperature (K), and mix.

        Includes evolving gases and particles, but no explicit surface point or
        short-lived species. Gas values are VMRs; particle values are condensed
        molecule equivalents per gas molecule. Raw values are not clipped.
        """
        sol = self.pc.mole_fraction_dict()
        pressure = sol['pressure']
        temperature = sol['temp']
        mix = {}
        species_names = self.pc.dat.species_names[:-2-self.pc.dat.nsl]
        for i,sp in enumerate(species_names):
            mix[sp] = sol[sp]
        return pressure, temperature, mix

    def initial_guess(self, edd, edd_pressure=None, surface_partial_pressures=None, **kwargs):

        # Extract surface pressure boundary conditions from photochemical model
        if surface_partial_pressures is None:
            pass

        # Create array of surface pressures for climate model
        P_i = np.ones(len(self.c.species_names))*1e-15
        for i,sp in enumerate(self.c.species_names):
            if sp in surface_partial_pressures:
                P_i[i] = surface_partial_pressures[sp]

        # Run the climate model and confirm convergence
        converged = self.c.RCE_robust(P_i)
        if not converged:
            raise Exception()

        # Extract results from climate run
        pressure, temperature, mix = self.profiles_from_climate()

        # Initialize photochemical code with climate results
        self.pc.initialize_atmosphere_p(
            pressure=pressure,
            temperature=temperature,
            edd=edd,
            mix=mix,
            particle_radius=None,
            persistent=True,
            maintain_toa_pressure=True,
            **kwargs
        )

        # Retrun
        return self.fixed_point_state_from_climate()
    
    def make_pressure_grid(self, log10P_surf):
        """Return log10 pressures from surface to climate top, in dyn/cm².

        Use one equally spaced log-pressure interval per climate layer,
        independent of the photochemical layer count.
        """
        top_pressure = self.c.P_top
        if not np.isfinite(top_pressure) or top_pressure <= 0.0:
            raise ValueError('Climate top pressure must be finite and positive (dyn/cm²).')
        log10P_top = np.log10(top_pressure)
        if not np.isfinite(log10P_surf) or log10P_surf <= log10P_top:
            raise ValueError('Surface pressure must be finite and exceed climate top pressure.')
        return np.linspace(log10P_surf, log10P_top, self.c.T.size + 1)

    def fixed_point_state_from_climate(self):
        """Return the fixed-point vector represented by the current climate.

        Neither model is modified, and climate convergence is not checked.
        Temperature is interpolated linearly in log10 pressure. Above the
        highest layer center, its temperature is held constant to ``P_top``.

        Returns
        -------
        numpy.ndarray
            ``[log10(P_surf), T_surface, ..., T_top]``, with pressure in
            dyn/cm² and temperatures in K. The temperature samples span
            equally spaced log10 pressures from the surface to climate top,
            with one interval per climate layer. The returned array owns its
            data and is independent of either model.
        """
        c = self.c
        pressure = np.concatenate(([c.P_surf], c.P))
        temperature = np.concatenate(([c.T_surf], c.T))
        if (pressure.ndim != 1 or temperature.shape != pressure.shape or
                pressure.size < 2):
            raise ValueError('Climate pressure and temperature profiles must have matching nonempty layers.')
        if not np.all(np.isfinite(pressure)) or np.any(pressure <= 0.0):
            raise ValueError('Climate pressures must be finite and positive; construct a climate profile first.')
        if np.any(np.diff(pressure) >= 0.0):
            raise ValueError('Climate pressure must strictly decrease from surface to top.')
        if not np.all(np.isfinite(temperature)) or np.any(temperature <= 0.0):
            raise ValueError('Climate temperatures must be finite and positive.')

        log10P_surf = np.log10(pressure[0])
        log10P = self.make_pressure_grid(log10P_surf)
        temperature_mapped = np.interp(
            log10P, np.log10(pressure[::-1]), temperature[::-1]
        )
        return np.concatenate(([log10P_surf], temperature_mapped))

    def g_eval(self, x, Pi, Kzz):
        c = self.c
        pc = self.pc

        log10P_surf, T = x[0], x[1:]
    
        #~~ Run photochemical model ~~#

        log10P = self.make_pressure_grid(log10P_surf)
        Kzz_ = np.ones(len(T))*Kzz
        pc.change_PT(10.0**log10P, T, Kzz_)
        # Integrate
        converged = pc.find_steady_state_robust()
        self.set_PTmix_from_photochem()

        if not converged:
            return np.ones_like(x)*np.nan

        # Setup climate model
        sol = pc.mole_fraction_dict()
        custom_dry_mix = {'pressure': sol['pressure']}
        P_i = np.ones(len(c.species_names))*1e-15
        for i,sp in enumerate(c.species_names):
            if sp in self.climate_condensates:
                P_i[c.species_names.index(sp)] = Pi[sp]
            else:
                P_i[c.species_names.index(sp)] = np.maximum(sol[sp][0],1e-30)*pc.var.surface_pressure*1e6
                custom_dry_mix[sp] = np.maximum(sol[sp],1e-30)
        
        # Solve climate (use last guess first)
        converged = c.RCE_robust(P_i, custom_dry_mix, T_init=np.append(c.T_surf, c.T), convecting_with_below_init=c.convecting_with_below)
        if not converged:
            converged = c.RCE_robust(P_i, custom_dry_mix)
        self.set_PTmix_from_climate()

        # Eject if not converged
        if not converged:
            return np.ones_like(x)*np.nan

        # Return
        result = self.fixed_point_state_from_climate()
        return result


class AdiabatClimateForCoupled(AdiabatClimate):

    def __init__(self, species_file, settings_file, flux_file, data_dir=None, top_pressure=1.0):

        # Climate model-top pressure in dyn/cm².
        top_pressure = float(top_pressure)
        if not np.isfinite(top_pressure) or top_pressure <= 0.0:
            raise ValueError('Climate top pressure must be finite and positive (dyn/cm²).')

        super().__init__(
            species_file, 
            settings_file, 
            flux_file,
            data_dir
        )

        self.solve_for_T_trop = True # Enable solving for T_trop
        self.max_rc_iters = 30 # Lots of iterations
        self.P_top = top_pressure
        self.convective_max_boundary_shift = 1
        self.max_rc_iters_convection = -1

    def surface_temperature_robust(self, P_i):

        T_guess_mid = self.rad.skin_temperature(0)*(1/2)**(-1/4) + 50
        T_perturbs = np.array([0.0, 20.0, -20.0, 30.0, 50.0, 80.0, 100.0, 200.0, 300.0, 400.0, 600.0])

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

    def RCE_robust(self, P_i, custom_dry_mix=None, T_init=None, convecting_with_below_init=None):
        
        if T_init is None:
            converged_simple = self.surface_temperature_robust(P_i)
            if not converged_simple:
                # If this fails, then we give up, returning no convergence
                return False
            T_surf_guess, T_guess, convecting_with_below_init = self.T_surf, self.T, self.convecting_with_below
        else:
            T_surf_guess = T_init[0]
            T_guess = T_init[1:].copy()

        try:
            converged = self.RCE(P_i, T_surf_guess, T_guess, convecting_with_below_init, custom_dry_mix)
        except ClimaException:
            converged = False

        return converged
