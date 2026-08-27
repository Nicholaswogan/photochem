import numpy as np
import numba as nb
from numba import types
from scipy import constants as const
from scipy import integrate
from tempfile import NamedTemporaryFile
import copy

from .._photochem import EvoAtmosphere
from .. import equilibrate
from ..utils._format import yaml, FormatSettings_main, MyDumper

###
### Extension of EvoAtmosphere class for gas giants
###

class GasGiantData():
    """Configuration and saved climate profiles for a gas-giant model.

    An instance is created as ``EvoAtmosphereGasGiant.gdat``. Pressure is
    in dyn/cm^2, temperature in K, eddy diffusion in cm^2/s, radius in cm, and
    mass in g. Attributes ending in ``_clima_grid`` retain the user's climate
    grid; ``P_desired``, ``T_desired``, and ``Kzz_desired`` contain the profile
    extended or truncated for photochemical initialization.
    """
    
    def __init__(self, planet_radius, planet_mass, P_ref, thermo_file):
        
        # Save several inputs
        self.planet_radius = planet_radius
        self.planet_mass = planet_mass
        self.P_ref = P_ref

        # Equilibrium solver
        self.gas = equilibrate.ChemEquiAnalysis(thermo_file)

        # Parameters using during initialization
        # The factor of pressure the atmosphere extends
        # compared to predicted quench points of gases
        self.BOA_pressure_factor = 5.0

        # If True, then the guessed initial condition will used
        # quenching relations as an initial guess
        self.initial_cond_with_quenching = True

        # Progress-reporting parameters.
        self.verbose = True # print information or not?
        self.freq_print = 100 # Frequency in which to print

        # Values that will be needed later. All of these set
        # in `initialize_to_climate_equilibrium_PT`
        self.P_clima_grid = None # The climate grid
        self.T_clima_grid = None # Prescribed temperature on the climate grid
        self.Kzz_clima_grid = None # Prescribed Kzz on the climate grid
        self.metallicity = None
        self.CtoO = None
        self.P_desired = None
        self.T_desired = None
        self.Kzz_desired = None
        # Index of climate grid that is bottom of photochemical grid
        self.ind_b = None

class EvoAtmosphereGasGiant(EvoAtmosphere):
    """Photochemical workflow for modeling gas-rich planets.

    The constructor configures static photochemical data. Call
    [initialize_to_climate_equilibrium_PT][photochem.extensions.EvoAtmosphereGasGiant.initialize_to_climate_equilibrium_PT] before integrating, or restore
    a previously saved state with [initialize_from_dict][photochem.extensions.EvoAtmosphereGasGiant.initialize_from_dict].

    Attributes
    ----------
    gdat : GasGiantData
        Gas-giant configuration and saved climate-grid profiles.
    dat, var, wrk
        Live views inherited from [photochem.EvoAtmosphere][].
    """

    def __init__(self, mechanism_file, stellar_flux_file, planet_mass, planet_radius, 
                 nz=100, photon_scale_factor=1.0, solar_zenith_angle=60.0, P_ref=1.0e6, 
                 thermo_file=None, data_dir=None):
        """Configure a gas-giant photochemical model.

        Parameters
        ----------
        mechanism_file : str
            Path to the file describing the reaction mechanism
        stellar_flux_file : str
            Path to the file describing the stellar UV flux.
        planet_mass : float
            Planet mass in grams
        planet_radius : float
            Planet radius in cm
        nz : int, optional
            The number of layers in the photochemical model, by default 100
        photon_scale_factor : float, optional
            Factor multiplying the input stellar photon flux, by default 1.
        solar_zenith_angle : float, optional
            Solar zenith angle in degrees, by default 60.
        P_ref : float, optional
            Pressure corresponding to ``planet_radius``, by default 1e6 dyn/cm^2.
        thermo_file : str, optional
            Optionally include a dedicated thermodynamic file.
        data_dir : str, optional
            Path to the data directory containing photolysis cross sections and other data
            needed to run the model. The packaged data are used by default.

        Notes
        -----
        Construction does not initialize an atmosphere. The pressure,
        temperature, eddy-diffusion, and composition profiles are supplied by
        a subsequent initialization method.
        """        
        
        # Configure the photochemical model. The atmosphere is initialized
        # later from the supplied climate profile.
        sol = yaml.safe_load(SETTINGS_TEMPLATE)
        sol['atmosphere-grid']['number-of-layers'] = int(nz)
        sol['planet']['planet-mass'] = float(planet_mass)
        sol['planet']['planet-radius'] = float(planet_radius)
        sol['planet']['photon-scale-factor'] = float(photon_scale_factor)
        sol['planet']['solar-zenith-angle'] = float(solar_zenith_angle)
        sol = FormatSettings_main(sol)
        with NamedTemporaryFile('w',suffix='.yaml') as f:
            yaml.dump(sol,f,Dumper=MyDumper)
            f.flush()
            super().__init__(
                mechanism_file,
                f.name,
                stellar_flux_file,
                data_dir=data_dir
            )

        if thermo_file is None:
            thermo_file = mechanism_file

        # New data
        self.gdat = GasGiantData(planet_radius, planet_mass, P_ref, thermo_file)

        # Values in photochem to adjust
        self.var.verbose = 0
        self.var.upwind_molec_diff = True
        self.var.atol = 1.0e-18
        self.var.conv_min_mix = 1e-10 # Min mix to consider during convergence check
        self.var.conv_longdy = 0.01 # threshold relative change that determines convergence
        self.var.custom_binary_diffusion_fcn = custom_binary_diffusion_fcn

        # Gas-giant defaults for the shared robust integration policy. These
        # preserve the extension's established limits while leaving one source
        # of truth for stepping, recovery, convergence, and TOA maintenance.
        self.var.nerrors_before_giveup = 10
        self.var.nsteps_before_conv_check = 300
        self.var.nsteps_before_reinit = 1000
        self.var.nsteps_before_giveup = 100_000
        maintenance = self.var.toa_pressure_maintenance
        maintenance.target_pressure = 0.1 # dyn/cm^2
        maintenance.pressure_factor = 3.0
        maintenance.nsteps_between_updates = 1000
        maintenance.max_failures = 2

    def initialize_to_climate_equilibrium_PT(self, P_in, T_in, Kzz_in, metallicity, CtoO, rainout_condensed_atoms=True):
        """Initialize from a climate profile and equilibrium composition.

        Input arrays must be one-dimensional, equal in length, and ordered
        from highest pressure (deepest) to lowest pressure. The initialized
        photochemical domain begins below the estimated quench levels and ends
        at ``var.toa_pressure_maintenance.target_pressure``.

        Parameters
        ----------
        P_in : ndarray, shape (nprofile,)
            Climate-grid pressure in dyn/cm^2, decreasing with index.
        T_in : ndarray, shape (nprofile,)
            Climate-grid temperature in K.
        Kzz_in : ndarray, shape (nprofile,)
            Climate-grid eddy diffusion in cm^2/s.
        metallicity : float
            Metallicity relative to solar.
        CtoO : float
            C/O ratio relative to solar. So CtoO = 1 is solar C/O ratio.
            CtoO = 2 is twice the solar C/O ratio.
        rainout_condensed_atoms : bool, optional
            Remove condensed atoms during equilibrium calculations, by default
            True.
        """

        gdat = self.gdat
        target_pressure = self._toa_pressure_target()

        if P_in.shape[0] != T_in.shape[0]:
            raise ValueError('P_in and T_in must have the same length.')
        if P_in.shape[0] != Kzz_in.shape[0]:
            raise ValueError('P_in and Kzz_in must have the same length.')

        # Save inputs
        gdat.P_clima_grid = P_in.copy()
        gdat.T_clima_grid = T_in.copy()
        gdat.Kzz_clima_grid = Kzz_in.copy()
        gdat.metallicity = metallicity
        gdat.CtoO = CtoO

        # Compute chemical equilibrium along the whole P-T profile
        mix, mubar = composition_at_metallicity(gdat.gas, T_in, P_in, CtoO, metallicity, rainout_condensed_atoms)

        # Altitude of P-T grid
        P1, T1, mubar1, z1 = compute_altitude_of_PT(P_in, gdat.P_ref, T_in, mubar, gdat.planet_radius, gdat.planet_mass, target_pressure)
        # If needed, extrapolate Kzz and mixing ratios
        if P1.shape[0] != Kzz_in.shape[0]:
            Kzz1 = np.append(Kzz_in,Kzz_in[-1])
            mix1 = {}
            for sp in mix:
                mix1[sp] = np.append(mix[sp],mix[sp][-1])
        else:
            Kzz1 = Kzz_in.copy()
            mix1 = mix

        # The gravity
        grav1 = gravity(gdat.planet_radius, gdat.planet_mass, z1)

        # Next, we compute the quench levels
        quench_levels = determine_quench_levels(T1, P1, Kzz1, mubar1, grav1)
        ind = np.min(quench_levels) # the deepest quench level

        # If desired, this bit applies quenched initial conditions, and recomputes
        # the altitude profile for this new mubar.
        if gdat.initial_cond_with_quenching:

            # Apply quenching to mixing ratios
            if "CH4" in mix1:
                mix1['CH4'][quench_levels[0]:] = mix1['CH4'][quench_levels[0]]
            if "CO" in mix1:
                mix1['CO'][quench_levels[0]:] = mix1['CO'][quench_levels[0]]
            if "NH3" in mix1:
                mix1['NH3'][quench_levels[2]:] = mix1['NH3'][quench_levels[2]]
            if "HCN" in mix1:
                mix1['HCN'][quench_levels[3]:] = mix1['HCN'][quench_levels[3]]
            if "H2" in mix1:
                # Quenching out H2 at the CH4 level seems to work well
                mix1['H2'][quench_levels[0]:] = mix1['H2'][quench_levels[0]]

            if "CO2" in mix1:
                # First, I need to equilibrate CO2 against quenched CO, H2O and H2.
                mix1['CO2'] = equilibrate_CO2_to_CO(mix1['CO'], mix1['H2O'], mix1['H2'], T1)
                # Next, I apply the quench.
                mix1['CO2'][quench_levels[1]:] = mix1['CO2'][quench_levels[1]]

            # Normalize mixing ratios
            mix_tot = np.zeros(mix1['H2'].shape[0])
            for key in mix1:
                mix_tot += mix1[key]
            for key in mix1:
                mix1[key] = mix1[key]/mix_tot

            # Compute mubar again
            mubar1[:] = 0.0
            for i,sp in enumerate(self.dat.species_names[:-2]):
                if sp in mix1:
                    for j in range(P1.shape[0]):
                        mubar1[j] += mix1[sp][j]*self.dat.species_mass[i]

            # Update z1 to get a new altitude profile
            P1, T1, mubar1, z1 = compute_altitude_of_PT(P1, gdat.P_ref, T1, mubar1, gdat.planet_radius, gdat.planet_mass, target_pressure)

        # Save the prescribed P-T-Kzz profile.
        gdat.P_desired = P1.copy()
        gdat.T_desired = T1.copy()
        gdat.Kzz_desired = Kzz1.copy()

        # Bottom of photochemical model will be at a pressure a factor
        # larger than the predicted quench pressure.
        if P1[ind]*gdat.BOA_pressure_factor > P1[0]:
            raise ValueError(
                'The climate pressure grid does not extend deeply enough for '
                'the requested photochemical lower boundary.'
            )
        gdat.ind_b = np.argmin(np.abs(P1 - P1[ind]*gdat.BOA_pressure_factor))
        
        self._initialize_atmosphere(P1, T1, Kzz1, z1, mix1)

    def reinitialize_to_new_climate_PT(self, P_in, T_in, Kzz_in, mix):
        """Reinitialize from updated climate profiles and composition.

        This method requires a prior call to
        [initialize_to_climate_equilibrium_PT][photochem.extensions.EvoAtmosphereGasGiant.initialize_to_climate_equilibrium_PT], and ``P_in`` must match
        that original climate pressure grid.

        Parameters
        ----------
        P_in : ndarray, shape (nprofile,)
            Climate-grid pressure in dyn/cm^2, decreasing with index.
        T_in : ndarray, shape (nprofile,)
            Climate-grid temperature in K.
        Kzz_in : ndarray, shape (nprofile,)
            Climate-grid eddy diffusion in cm^2/s.
        mix : dict[str, ndarray]
            Mixing-ratio profile for every gas species. Each value has shape
            ``(nprofile,)``.

        """

        gdat = self.gdat
        target_pressure = self._toa_pressure_target()

        if gdat.P_clima_grid is None:
            raise RuntimeError(
                'reinitialize_to_new_climate_PT requires a prior call to '
                'initialize_to_climate_equilibrium_PT.'
            )
        if not np.all(np.isclose(gdat.P_clima_grid,P_in)):
            raise ValueError('P_in does not match the saved climate pressure grid.')
        if P_in.shape[0] != T_in.shape[0]:
            raise ValueError('P_in and T_in must have the same length.')
        if P_in.shape[0] != Kzz_in.shape[0]:
            raise ValueError('P_in and Kzz_in must have the same length.')
        for key in mix:
            if P_in.shape[0] != mix[key].shape[0]:
                raise ValueError(
                    f'The mixing-ratio profile for {key!r} must have the '
                    'same length as P_in.'
                )
        # Require all gases be specified. Particles can be ignored.
        if set(list(mix.keys())) != set(self.dat.species_names[self.dat.np:(-2-self.dat.nsl)]):
            raise ValueError('mix must contain every gas species and no extra species.')

        # Save the newly prescribed climate-grid profiles. The pressure grid
        # is required to match the one saved during initial initialization.
        gdat.T_clima_grid = T_in.copy()
        gdat.Kzz_clima_grid = Kzz_in.copy()
        
        # Compute mubar
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        mubar = np.zeros(T_in.shape[0])
        species_mass = self.dat.species_mass
        particle_names = self.dat.species_names[:self.dat.np]
        for sp in mix:
            if sp not in particle_names:
                ind = species_names.index(sp)
                mubar = mubar + mix[sp]*species_mass[ind]

        # Compute altitude of P-T grid
        P1, T1, mubar1, z1 = compute_altitude_of_PT(P_in, gdat.P_ref, T_in, mubar, gdat.planet_radius, gdat.planet_mass, target_pressure)
        # If needed, extrapolte Kzz and mixing ratios
        if P1.shape[0] != Kzz_in.shape[0]:
            Kzz1 = np.append(Kzz_in,Kzz_in[-1])
            mix1 = {}
            for sp in mix:
                mix1[sp] = np.append(mix[sp],mix[sp][-1])
        else:
            Kzz1 = Kzz_in.copy()
            mix1 = mix

        # Save the prescribed P-T-Kzz profile.
        gdat.P_desired = P1.copy()
        gdat.T_desired = T1.copy()
        gdat.Kzz_desired = Kzz1.copy()

        self._initialize_atmosphere(P1, T1, Kzz1, z1, mix1)

    def _initialize_atmosphere(self, P1, T1, Kzz1, z1, mix1):
        "Initialize the shared photochemical model from gas-giant profiles."

        gdat = self.gdat
        target_pressure = self._toa_pressure_target()

        # Select the pressure interval used by the photochemical model. If the
        # climate grid continues above the requested model top, interpolate an
        # exact endpoint rather than selecting a climate level on an arbitrary
        # side of the target. The complete climate profile remains in gdat.
        P1 = np.asarray(P1)
        T1 = np.asarray(T1)
        Kzz1 = np.asarray(Kzz1)
        z1 = np.asarray(z1)
        top_candidates = np.flatnonzero(P1 <= target_pressure)
        if top_candidates.size == 0:
            raise ValueError(
                'The supplied pressure profile does not reach the requested '
                'photochemical model top.'
            )
        ind_t = top_candidates[0]
        if ind_t <= gdat.ind_b:
            raise ValueError('The requested photochemical pressure domain is empty.')

        if np.isclose(P1[ind_t], target_pressure, rtol=1.0e-12, atol=0.0):
            inds = slice(gdat.ind_b, ind_t + 1)
            pressure_profile = P1[inds]
            temperature_profile = T1[inds]
            edd_profile = Kzz1[inds]
            mix_profile = {
                sp: np.asarray(values)[inds]
                for sp, values in mix1.items()
            }
        else:
            inds = slice(gdat.ind_b, ind_t)
            log_pressure = np.log(P1[::-1])
            log_target = np.log(target_pressure)
            pressure_profile = np.append(P1[inds], target_pressure)
            temperature_profile = np.append(
                T1[inds], np.interp(log_target, log_pressure, T1[::-1])
            )
            edd_profile = np.append(
                Kzz1[inds],
                10.0**np.interp(
                    log_target, log_pressure, np.log10(Kzz1[::-1])
                )
            )
            mix_profile = {}
            for sp, values in mix1.items():
                values = np.asarray(values)
                value_at_top = 10.0**np.interp(
                    log_target, log_pressure,
                    np.log10(np.clip(values[::-1], 1.0e-100, np.inf))
                )
                mix_profile[sp] = np.append(values[inds], value_at_top)

        # The base class interprets planet_radius at the lower pressure
        # boundary. Convert from the gas-giant reference pressure before
        # delegating the hydrostatic construction and profile mapping.
        planet_radius_new = gdat.planet_radius + z1[gdat.ind_b]
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        mix_profile = {
            sp: values for sp, values in mix_profile.items()
            if sp in species_names
        }

        planet_radius_old = self.dat.planet_radius
        self.dat.planet_radius = planet_radius_new
        try:
            self.initialize_atmosphere_p(
                pressure_profile, temperature_profile, edd_profile,
                mix_profile, persistent=True,
                maintain_toa_pressure=True,
                target_pressure=target_pressure
            )
        except Exception:
            self.dat.planet_radius = planet_radius_old
            raise
        # Retain the gas-giant boundary-condition policy. Pressure boundary
        # conditions use the initialized lowest model cell, as before.
        for i,sp in enumerate(species_names):
            if i >= self.dat.np:
                self.set_lower_bc(sp, bc_type='Moses') # gas
            else:
                self.set_lower_bc(sp, bc_type='vdep', vdep=0.0) # particle
        particle_names = self.dat.species_names[:self.dat.np]
        for sp in mix_profile:
            if sp not in particle_names:
                ind = species_names.index(sp)
                mix_surf = self.wrk.usol[ind,0]/self.wrk.density[0]
                Pi = self.wrk.pressure_hydro[0]*mix_surf
                self.set_lower_bc(sp, bc_type='press', press=Pi)

    def _toa_pressure_target(self):
        "Return the shared, validated gas-giant TOA-pressure target."
        target_pressure = self.var.toa_pressure_maintenance.target_pressure
        if not np.isfinite(target_pressure) or target_pressure <= 0.0:
            raise ValueError(
                'var.toa_pressure_maintenance.target_pressure must be finite and positive'
            )
        return target_pressure

    def return_atmosphere_climate_grid(self):
        """Return prescribed profiles and composition on the climate grid.

        Temperature and eddy diffusion are the saved climate inputs. Species
        mixing ratios are log-interpolated or constant-extrapolated from the
        photochemical result.

        Returns
        -------
        dict[str, ndarray]
            ``pressure`` (dyn/cm^2), ``temperature`` (K), ``Kzz`` (cm^2/s),
            and a mixing-ratio array for each modeled species.
        """

        gdat = self.gdat

        if gdat.P_clima_grid is None:
            raise RuntimeError(
                'return_atmosphere_climate_grid requires an initialized '
                'gas-giant atmosphere.'
            )

        # return full atmosphere
        out = self.return_atmosphere()

        # Temperature and Kzz are prescribed climate-model inputs, so return
        # them directly. Only the photochemical composition needs mapping to
        # the climate grid.
        sol = {}
        sol['pressure'] = gdat.P_clima_grid.copy()
        sol['temperature'] = gdat.T_clima_grid.copy()
        sol['Kzz'] = gdat.Kzz_clima_grid.copy()
        log10Pclima = np.log10(gdat.P_clima_grid[::-1]).copy()
        log10P = np.log10(out['pressure'][::-1]).copy()

        for key in out:
            if key not in ['pressure','temperature','Kzz']:
                tmp = np.log10(np.clip(out[key][::-1].copy(),a_min=1e-100,a_max=np.inf))
                mix = np.interp(log10Pclima, log10P, tmp)
                sol[key] = 10.0**mix[::-1].copy()

        return sol

    def return_atmosphere(self, include_deep_atmosphere = True, equilibrium = False, rainout_condensed_atoms = True):
        """Return atmospheric profiles on the photochemical pressure grid.

        Parameters
        ----------
        include_deep_atmosphere : bool, optional
            Prepend the deeper prescribed atmosphere, evaluated in chemical
            equilibrium, by default True.
        equilibrium : bool, optional
            Return equilibrium rather than photochemical mixing ratios on the
            photochemical grid, by default False.
        rainout_condensed_atoms : bool, optional
            Remove condensed atoms in equilibrium calculations, by default
            True.

        Returns
        -------
        dict[str, ndarray]
            ``pressure`` (dyn/cm^2), ``temperature`` (K), ``Kzz`` (cm^2/s),
            and a mixing-ratio array for each modeled species.
        """

        gdat = self.gdat      

        if gdat.P_clima_grid is None:
            raise RuntimeError(
                'return_atmosphere requires an initialized gas-giant atmosphere.'
            )

        out = {}
        out['pressure'] = self.wrk.pressure_hydro
        out['temperature'] = self.var.temperature
        out['Kzz'] = self.var.edd
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        if equilibrium:
            mix, mubar = composition_at_metallicity(gdat.gas, out['temperature'], out['pressure'], gdat.CtoO, gdat.metallicity, rainout_condensed_atoms)
            for key in mix:
                out[key] = mix[key]
            for key in species_names[:self.dat.np]:
                out[key] = np.zeros(mix['H2'].shape[0])
        else:
            for i,sp in enumerate(species_names):
                mix = self.wrk.usol[i,:]/self.wrk.density
                out[sp] = mix

        if not include_deep_atmosphere:
            return out

        # Prepend the deeper atmosphere, which we will assume is at Equilibrium
        inds = np.where(gdat.P_desired > self.wrk.pressure_hydro[0])
        out1 = {}
        out1['pressure'] = gdat.P_desired[inds]
        out1['temperature'] = gdat.T_desired[inds]
        out1['Kzz'] = gdat.Kzz_desired[inds]
        mix, mubar = composition_at_metallicity(gdat.gas, out1['temperature'], out1['pressure'], gdat.CtoO, gdat.metallicity, rainout_condensed_atoms)
        
        out['pressure'] = np.append(out1['pressure'],out['pressure'])
        out['temperature'] = np.append(out1['temperature'],out['temperature'])
        out['Kzz'] = np.append(out1['Kzz'],out['Kzz'])
        for i,sp in enumerate(species_names):
            if sp in mix:
                out[sp] = np.append(mix[sp],out[sp])
            else:
                out[sp] = np.append(np.zeros(mix['H2'].shape[0]),out[sp])

        return out
    
    def robust_step(self):
        """Take one shared robust step and optionally report gas-giant progress.

        Stepping, recovery, convergence, restart, limit, and TOA-maintenance
        policy are implemented by [EvoAtmosphere][photochem.EvoAtmosphere]. This override only
        retains the extension's compact progress display.

        Returns
        -------
        tuple
            Two booleans ``give_up, reached_steady_state`` from the shared
            robust stepper.
        """
        nsteps_before = self.wrk.nsteps_total
        give_up, reached_steady_state = super().robust_step()

        nsteps = self.wrk.nsteps_total
        accepted_step = nsteps > nsteps_before
        report_interval = max(int(self.gdat.freq_print), 1)
        if self.gdat.verbose and (
            reached_steady_state or give_up or
            (accepted_step and not (nsteps % report_interval))
        ):
            TOA_pressure = self.wrk.pressure_hydro[-1]
            print('nsteps = %i  longdy = %.1e  TOA_pressure = %.1e' %
                  (nsteps, self.wrk.longdy, TOA_pressure/1e6))

        return give_up, reached_steady_state
    
    def model_state_to_dict(self):
        """Serialize the initialized gas-giant atmospheric state.

        Returns
        -------
        dict
            State accepted by [initialize_from_dict][photochem.extensions.EvoAtmosphereGasGiant.initialize_from_dict].
        """

        gdat = self.gdat

        if gdat.P_clima_grid is None:
            raise RuntimeError(
                'model_state_to_dict requires an initialized gas-giant atmosphere.'
            )

        out = {}
        out['P_clima_grid'] = gdat.P_clima_grid
        out['T_clima_grid'] = gdat.T_clima_grid
        out['Kzz_clima_grid'] = gdat.Kzz_clima_grid
        out['metallicity'] = gdat.metallicity
        out['CtoO'] = gdat.CtoO
        out['P_desired'] = gdat.P_desired
        out['T_desired'] = gdat.T_desired
        out['Kzz_desired'] = gdat.Kzz_desired
        out['toa_pressure_target'] = self._toa_pressure_target()
        out['ind_b'] = gdat.ind_b
        out['planet_radius_new'] = self.dat.planet_radius
        out['top_atmos'] = self.var.top_atmos
        out['temperature'] = self.var.temperature
        out['edd'] = self.var.edd
        out['usol'] = self.wrk.usol
        out['P_i_surf'] = (self.wrk.usol[self.dat.np:,0]/self.wrk.density[0])*self.wrk.pressure[0]

        return out

    def initialize_from_dict(self, out):
        """Restore a state produced by [model_state_to_dict][photochem.extensions.EvoAtmosphereGasGiant.model_state_to_dict].

        Any active stepper and pressure-profile prescription are replaced.

        Parameters
        ----------
        out : dict
            Serialized gas-giant state.
        """

        gdat = self.gdat
        target_pressure = out.get('toa_pressure_target', 0.1)
        if not np.isfinite(target_pressure) or target_pressure <= 0.0:
            raise ValueError('Saved TOA-pressure target must be finite and positive')

        # The saved state replaces any current integration and prescribed
        # pressure profile.
        self.destroy_stepper()
        self.clear_press_temp_edd_profile()

        gdat.P_clima_grid = out['P_clima_grid']
        # Fall back to the leading prescribed-profile values for states saved
        # before the climate-grid temperature and Kzz fields were introduced.
        nclima = gdat.P_clima_grid.shape[0]
        gdat.T_clima_grid = out.get(
            'T_clima_grid', out['T_desired'][:nclima]
        )
        gdat.Kzz_clima_grid = out.get(
            'Kzz_clima_grid', out['Kzz_desired'][:nclima]
        )
        gdat.metallicity = out['metallicity']
        gdat.CtoO = out['CtoO']
        gdat.P_desired = out['P_desired']
        gdat.T_desired = out['T_desired']
        gdat.Kzz_desired = out['Kzz_desired']
        gdat.ind_b = out['ind_b']
        self.dat.planet_radius = out['planet_radius_new']
        self.update_vertical_grid(TOA_alt=out['top_atmos'])
        self.set_temperature(out['temperature'])
        self.var.edd = out['edd']
        self.wrk.usol = out['usol']

        # Now set boundary conditions
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        for i,sp in enumerate(species_names):
            if i >= self.dat.np:
                self.set_lower_bc(sp, bc_type='Moses') # gas
            else:
                self.set_lower_bc(sp, bc_type='vdep', vdep=0.0) # particle
        species_names = self.dat.species_names[self.dat.np:(-2-self.dat.nsl)]
        for i,sp in enumerate(species_names):
            self.set_lower_bc(sp, bc_type='press', press=out['P_i_surf'][i])

        self.set_press_temp_edd_profile(
            gdat.P_desired, gdat.T_desired, gdat.Kzz_desired,
            hydro_pressure=True, maintain_toa_pressure=True,
            target_pressure=target_pressure
        )

###
### Helper functions for the EvoAtmosphereGasGiant class
###

@nb.cfunc(nb.double(nb.double, nb.double, nb.double))
def custom_binary_diffusion_fcn(mu_i, mubar, T):
    # Equation 6 in Gladstone et al. (1996)
    b = 3.64e-5*T**(1.75-1.0)*7.3439e21*np.sqrt(2.01594/mu_i)
    return b

@nb.njit()
def CH4_CO_quench_timescale(T, P):
    "T in K, P in dyn/cm^2, tq in s. Equation 11."
    P_bars = P/1.0e6
    tq = 3.0e-6*P_bars**-1*np.exp(42_000.0/T)
    return tq

@nb.njit()
def NH3_quench_timescale(T, P):
    "T in K, P in dyn/cm^2, tq in s. Equation 32."
    P_bars = P/1.0e6
    tq = 1.0e-7*P_bars**-1*np.exp(52_000.0/T)
    return tq

@nb.njit()
def HCN_quench_timescale(T, P):
    "T in K, P in dyn/cm^2, tq in s. From PICASO."
    P_bars = P/1.0e6
    tq = (1.5e-4/(P_bars*(3.0**0.7)))*np.exp(36_000.0/T)
    return tq

@nb.njit()
def CO2_quench_timescale(T, P):
    "T in K, P in dyn/cm^2, tq in s. Equation 44."
    P_bars = P/1.0e6
    tq = 1.0e-10*P_bars**-0.5*np.exp(38_000.0/T)
    return tq

@nb.njit()
def equilibrate_CO2_to_CO(fCO, fH2O, fH2, T):
    """The mole fraction of CO2 in equilibrium with CO, H2O and H2.
    From Equation 43 in Zahnle and Marley (2014).

    Parameters
    ----------
    fCO : float
        The CO mole fraction
    fH2O : float
        The H2O mole fraction
    fH2 : float
        The H2 mole fraction
    T : float
        Temperature in K

    Returns
    -------
    float
        The CO2 mole fraction
    """    
    K = 18.3*np.exp(-2376/T - (932/T)**2)
    fCO2 = (fCO*fH2O)/(K*fH2)
    return fCO2

@nb.njit()
def scale_height(T, mubar, grav):
    "All inputs are CGS."
    k_boltz = const.k*1e7
    H = (const.Avogadro*k_boltz*T)/(mubar*grav)
    return H

@nb.njit()
def determine_quench_levels(T, P, Kzz, mubar, grav):

    # Mixing timescale
    tau_mix = scale_height(T, mubar, grav)**2/Kzz

    # Quenching timescales
    tau_CH4 = CH4_CO_quench_timescale(T, P)
    tau_CO2 = CO2_quench_timescale(T, P)
    tau_NH3 = NH3_quench_timescale(T, P)
    tau_HCN = HCN_quench_timescale(T, P)

    # Quench level is when the chemistry timescale
    # exceeds the mixing timescale.
    quench_levels = np.zeros(4, dtype=np.int32)
    
    for i in range(P.shape[0]):
        quench_levels[0] = i
        if tau_CH4[i] > tau_mix[i]:
            break

    for i in range(P.shape[0]):
        quench_levels[1] = i
        if tau_CO2[i] > tau_mix[i]:
            break

    for i in range(P.shape[0]):
        quench_levels[2] = i
        if tau_NH3[i] > tau_mix[i]:
            break

    for i in range(P.shape[0]):
        quench_levels[3] = i
        if tau_HCN[i] > tau_mix[i]:
            break

    return quench_levels
    
@nb.experimental.jitclass()
class TempPressMubar:

    log10P : types.double[:] # type: ignore
    T : types.double[:] # type: ignore
    mubar : types.double[:] # type: ignore

    def __init__(self, P, T, mubar):
        self.log10P = np.log10(P)[::-1].copy()
        self.T = T[::-1].copy()
        self.mubar = mubar[::-1].copy()

    def temperature_mubar(self, P):
        T = np.interp(np.log10(P), self.log10P, self.T)
        mubar = np.interp(np.log10(P), self.log10P, self.mubar)
        return T, mubar

@nb.njit()
def gravity(radius, mass, z):
    G_grav = const.G
    grav = G_grav * (mass/1.0e3) / ((radius + z)/1.0e2)**2.0
    grav = grav*1.0e2 # convert to cgs
    return grav

@nb.njit()
def hydrostatic_equation(P, u, planet_radius, planet_mass, ptm):
    z = u[0]
    grav = gravity(planet_radius, planet_mass, z)
    T, mubar = ptm.temperature_mubar(P)
    k_boltz = const.Boltzmann*1e7
    dz_dP = -(k_boltz*T*const.Avogadro)/(mubar*grav*P)
    return np.array([dz_dP])

def compute_altitude_of_PT(P, P_ref, T, mubar, planet_radius, planet_mass, P_top):
    ptm = TempPressMubar(P, T, mubar)
    args = (planet_radius, planet_mass, ptm)

    if P_top < P[-1]:
        # If P_top is lower P than P grid, then we extend it
        P_top_ = P_top
        P_ = np.append(P,P_top_)
        T_ = np.append(T,T[-1])
        mubar_ = np.append(mubar,mubar[-1])
    else:
        P_top_ = P[-1]
        P_ = P.copy()
        T_ = T.copy()
        mubar_ = mubar.copy()

    # Make sure P_ref is in the P grid
    if P_ref > P_[0] or P_ref < P_[-1]:
        raise ValueError('P_ref must lie within the pressure grid.')
    
    # Find first index with lower pressure than P_ref
    ind = 0
    for i in range(P_.shape[0]):
        if P_[i] < P_ref:
            ind = i
            break

    # Integrate from P_ref to TOA
    out2 = integrate.solve_ivp(hydrostatic_equation, [P_ref, P_[-1]], np.array([0.0]), t_eval=P_[ind:], args=args, rtol=1e-6)
    if not out2.success:
        raise RuntimeError(
            'Hydrostatic integration from P_ref to the model top failed: '
            + out2.message
        )
    # Integrate from P_ref to BOA
    out1 = integrate.solve_ivp(hydrostatic_equation, [P_ref, P_[0]], np.array([0.0]), t_eval=P_[:ind][::-1], args=args, rtol=1e-6)
    if not out1.success:
        raise RuntimeError(
            'Hydrostatic integration from P_ref to the model bottom failed: '
            + out1.message
        )

    # Stitch together
    z_ = np.append(out1.y[0][::-1],out2.y[0])

    return P_, T_, mubar_, z_

###
### A simple metallicity calculator
###

def composition_at_metallicity(gas, T, P, CtoO, metal, rainout_condensed_atoms = True):
    """Given a T-P profile, C/O ratio and metallicity, the code
    computes chemical equilibrium composition.

    Parameters
    ----------
    gas : ChemEquiAnalysis
    T : ndarray[dim=1,float64]
        Temperature in K
    P : ndarray[dim=1,float64]
        Pressure in dyn/cm^2.
    CtoO : float
        The C / O ratio relative to solar. CtoO = 1 would be the same
        composition as solar.
    metal : float
        Metallicity relative to solar.
    rainout_condensed_atoms : bool, optional
        If True, then the code will rainout atoms that condense.

    Returns
    -------
    dict
        Composition at chemical equilibrium.
    """

    # Check T and P
    if isinstance(T, float) or isinstance(T, int):
        T = np.array([T],np.float64)
    if isinstance(P, float) or isinstance(P, int):
        P = np.array([P],np.float64)
    if not isinstance(P, np.ndarray):
        raise ValueError('P must be a NumPy array.')
    if not isinstance(T, np.ndarray):
        raise ValueError('T must be a NumPy array.')
    if T.ndim != 1:
        raise ValueError('"T" must have one dimension')
    if P.ndim != 1:
        raise ValueError('"P" must have one dimension')
    if T.shape[0] != P.shape[0]:
        raise ValueError('"P" and "T" must be the same length')
    # Check CtoO and metal
    if CtoO <= 0:
        raise ValueError('"CtoO" must be greater than 0')
    if metal <= 0:
        raise ValueError('"metal" must be greater than 0')

    # For output
    out = {}
    for sp in gas.gas_names:
        out[sp] = np.empty(P.shape[0])
    mubar = np.empty(P.shape[0])
    
    molfracs_atoms = gas.molfracs_atoms_sun
    for i,sp in enumerate(gas.atoms_names):
        if sp != 'H' and sp != 'He':
            molfracs_atoms[i] = gas.molfracs_atoms_sun[i]*metal
    molfracs_atoms = molfracs_atoms/np.sum(molfracs_atoms)

    # Adjust C and O to get desired C/O ratio. CtoO is relative to solar
    indC = gas.atoms_names.index('C')
    indO = gas.atoms_names.index('O')
    x = CtoO*(molfracs_atoms[indC]/molfracs_atoms[indO])
    a = (x*molfracs_atoms[indO] - molfracs_atoms[indC])/(1+x)
    molfracs_atoms[indC] = molfracs_atoms[indC] + a
    molfracs_atoms[indO] = molfracs_atoms[indO] - a

    # Compute chemical equilibrium at all altitudes
    for i in range(P.shape[0]):
        if i > 0:
            gas.use_prev_guess = True
        for eps in [0.0, 1.0e-12, -1.0e-12, 1.0e-8, -1.0e-8]:
            converged = gas.solve(P[i], T[i] + T[i]*eps, molfracs_atoms=molfracs_atoms)
            if converged:
                break
        # Do not enforce convergence.
        for j,sp in enumerate(gas.gas_names):
            out[sp][i] = gas.molfracs_species_gas[j]
        mubar[i] = gas.mubar
        if rainout_condensed_atoms:
            molfracs_atoms = gas.molfracs_atoms_gas
    gas.use_prev_guess = False

    return out, mubar

###
### Template input files for Photochem
###

SETTINGS_TEMPLATE = \
"""
atmosphere-grid:
  number-of-layers: NULL

planet:
  planet-mass: NULL
  planet-radius: NULL
  surface-albedo: 0.0
  solar-zenith-angle: 60.0
  hydrogen-escape:
    type: none
  default-gas-lower-boundary: Moses
  water:
    gas-rainout: false

boundary-conditions:
- name: H2
  lower-boundary: {type: Moses}
  upper-boundary: {type: veff, veff: 0}
"""
