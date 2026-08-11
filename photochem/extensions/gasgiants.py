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
        self.metallicity = None
        self.CtoO = None
        self.P_desired = None
        self.T_desired = None
        self.Kzz_desired = None
        # Index of climate grid that is bottom of photochemical grid
        self.ind_b = None

class EvoAtmosphereGasGiant(EvoAtmosphere):
    "An extension to the EvoAtmosphere class for modeling gas-rich planets."

    def __init__(self, mechanism_file, stellar_flux_file, planet_mass, planet_radius, 
                 nz=100, photon_scale_factor=1.0, solar_zenith_angle=60.0, P_ref=1.0e6, 
                 thermo_file=None, data_dir=None):
        """Initializes the code

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
            description
        solar_zenith_angle : float, optional
            description
        P_ref : float, optional
            Pressure level corresponding to the planet_radius, by default 1e6 dynes/cm^2
        thermo_file : str, optional
            Optionally include a dedicated thermodynamic file.
        data_dir : str, optional
            Path to the data directory containing photolysis cross sections and other data
            needed to run the model
        """        
        
        # Configure the photochemical model. The atmosphere is initialized
        # later from the supplied pressure-based climate profile.
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
        self.var.autodiff = True # Turn on autodiff
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
        maintenance.target_pressure = 0.1 # dynes/cm^2
        maintenance.pressure_factor = 3.0
        maintenance.nsteps_between_updates = 1000
        maintenance.max_failures = 2

    def initialize_to_climate_equilibrium_PT(self, P_in, T_in, Kzz_in, metallicity, CtoO, rainout_condensed_atoms=True):
        """Initialized the photochemical model to a climate model result that assumes chemical equilibrium
        at some metallicity and C/O ratio.

        Parameters
        ----------
        P_in : ndarray[dim=1,double]
            The pressures in the climate grid (dynes/cm^2). P_in[0] is pressure at
            the deepest layer of the atmosphere
        T_in : ndarray[dim1,double]
            The temperatures in the climate grid corresponding to P_in (K)
        Kzz_in : ndarray[dim1,double]
            The eddy diffusion at each pressure P_in (cm^2/s)
        metallicity : float
            Metallicity relative to solar.
        CtoO : float
            C/O ratio relative to solar. So CtoO = 1 is solar C/O ratio.
            CtoO = 2 is twice the solar C/O ratio.
        """

        gdat = self.gdat
        target_pressure = self._toa_pressure_target()

        if P_in.shape[0] != T_in.shape[0]:
            raise Exception('Input P and T must have same shape')
        if P_in.shape[0] != Kzz_in.shape[0]:
            raise Exception('Input P and Kzz must have same shape')

        # Save inputs
        gdat.P_clima_grid = P_in
        gdat.metallicity = metallicity
        gdat.CtoO = CtoO

        # Compute chemical equilibrium along the whole P-T profile
        mix, mubar = composition_at_metallicity(gdat.gas, T_in, P_in, CtoO, metallicity, rainout_condensed_atoms)

        if target_pressure*3 > P_in[-1]:
            raise Exception('The photochemical grid needs to extend above the climate grid')

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
            raise Exception('BOA in photochemical model wants to be deeper than BOA of climate model.')
        gdat.ind_b = np.argmin(np.abs(P1 - P1[ind]*gdat.BOA_pressure_factor))
        
        self._initialize_atmosphere(P1, T1, Kzz1, z1, mix1)

    def reinitialize_to_new_climate_PT(self, P_in, T_in, Kzz_in, mix):
        """Reinitializes the photochemical model to the input P, T, Kzz, and mixing ratios
        from the climate model.

        Parameters
        ----------
        P_in : ndarray[ndim=1,double]
            Pressure grid in climate model (dynes/cm^2).
        T_in : ndarray[ndim=1,double]
            Temperatures corresponding to P_in (K)
        Kzz_in : ndarray[ndim,double]
            Eddy diffusion coefficients at each pressure level (cm^2/s)
        mix : dict
            Mixing ratios of all species in the atmosphere

        """

        gdat = self.gdat
        target_pressure = self._toa_pressure_target()

        if gdat.P_clima_grid is None:
            raise Exception('This routine can only be called after `initialize_to_climate_equilibrium_PT`')
        if not np.all(np.isclose(gdat.P_clima_grid,P_in)):
            raise Exception('Input pressure grid does not match saved pressure grid')
        if P_in.shape[0] != T_in.shape[0]:
            raise Exception('Input P and T must have same shape')
        if P_in.shape[0] != Kzz_in.shape[0]:
            raise Exception('Input P and Kzz must have same shape')
        for key in mix:
            if P_in.shape[0] != mix[key].shape[0]:
                raise Exception('Input P and mix must have same shape')
        # Require all gases be specified. Particles can be ignored.
        if set(list(mix.keys())) != set(self.dat.species_names[self.dat.np:(-2-self.dat.nsl)]):
            raise Exception('Some species are missing from input mix') 
        
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

        # Select the pressure interval used by the photochemical model. The
        # deeper profile remains in gdat for return_atmosphere.
        ind_t = np.argmin(np.abs(P1 - target_pressure))
        if ind_t <= gdat.ind_b:
            raise Exception('The photochemical pressure domain is empty.')
        inds = slice(gdat.ind_b, ind_t + 1)

        # The base class interprets planet_radius at the lower pressure
        # boundary. Convert from the gas-giant reference pressure before
        # delegating the hydrostatic construction and profile mapping.
        planet_radius_new = gdat.planet_radius + z1[gdat.ind_b]
        species_names = self.dat.species_names[:(-2-self.dat.nsl)]
        mix_profile = {
            sp: np.asarray(mix1[sp])[inds]
            for sp in mix1 if sp in species_names
        }

        planet_radius_old = self.dat.planet_radius
        self.dat.planet_radius = planet_radius_new
        try:
            self.initialize_atmosphere_p(
                np.asarray(P1)[inds], np.asarray(T1)[inds],
                np.asarray(Kzz1)[inds], mix_profile, persistent=True,
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
        """Returns a dictionary with temperature, Kzz and mixing ratios
        on the climate model grid.

        Returns
        -------
        dict
            Contains temperature, Kzz, and mixing ratios.
        """

        gdat = self.gdat

        if gdat.P_clima_grid is None:
            raise Exception('This routine can only be called after `initialize_to_climate_equilibrium_PT`')

        # return full atmosphere
        out = self.return_atmosphere()

        # Interpolate full atmosphere to clima grid
        sol = {}
        sol['pressure'] = gdat.P_clima_grid.copy()
        log10Pclima = np.log10(gdat.P_clima_grid[::-1]).copy()
        log10P = np.log10(out['pressure'][::-1]).copy()

        T = np.interp(log10Pclima, log10P, out['temperature'][::-1].copy())
        sol['temperature'] = T[::-1].copy()

        Kzz = np.interp(log10Pclima, log10P, np.log10(out['Kzz'][::-1].copy()))
        sol['Kzz'] = 10.0**Kzz[::-1].copy()

        for key in out:
            if key not in ['pressure','temperature','Kzz']:
                tmp = np.log10(np.clip(out[key][::-1].copy(),a_min=1e-100,a_max=np.inf))
                mix = np.interp(log10Pclima, log10P, tmp)
                sol[key] = 10.0**mix[::-1].copy()

        return sol

    def return_atmosphere(self, include_deep_atmosphere = True, equilibrium = False, rainout_condensed_atoms = True):
        """Returns a dictionary with temperature, Kzz and mixing ratios
        on the photochemical grid.

        Parameters
        ----------
        include_deep_atmosphere : bool, optional
            If True, then results will include portions of the deep
            atomsphere that are not part of the photochemical grid, by default True

        Returns
        -------
        dict
            Contains temperature, Kzz, and mixing ratios.
        """

        gdat = self.gdat      

        if gdat.P_clima_grid is None:
            raise Exception('This routine can only be called after `initialize_to_climate_equilibrium_PT`')

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
        policy are implemented by :class:`EvoAtmosphere`. This override only
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
        """Returns a dictionary containing all information needed to reinitialize the atmospheric
        state. This dictionary can be used as an input to "initialize_from_dict".
        """

        gdat = self.gdat

        if gdat.P_clima_grid is None:
            raise Exception('This routine can only be called after `initialize_to_climate_equilibrium_PT`')

        out = {}
        out['P_clima_grid'] = gdat.P_clima_grid
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
        """Initializes the model from a dictionary created by the "model_state_to_dict" routine.
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
    "T in K, P in dynes/cm^2, tq in s. Equation 11."
    P_bars = P/1.0e6
    tq = 3.0e-6*P_bars**-1*np.exp(42_000.0/T)
    return tq

@nb.njit()
def NH3_quench_timescale(T, P):
    "T in K, P in dynes/cm^2, tq in s. Equation 32."
    P_bars = P/1.0e6
    tq = 1.0e-7*P_bars**-1*np.exp(52_000.0/T)
    return tq

@nb.njit()
def HCN_quench_timescale(T, P):
    "T in K, P in dynes/cm^2, tq in s. From PICASO."
    P_bars = P/1.0e6
    tq = (1.5e-4/(P_bars*(3.0**0.7)))*np.exp(36_000.0/T)
    return tq

@nb.njit()
def CO2_quench_timescale(T, P):
    "T in K, P in dynes/cm^2, tq in s. Equation 44."
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
        raise Exception('Reference pressure must be within P grid.')
    
    # Find first index with lower pressure than P_ref
    ind = 0
    for i in range(P_.shape[0]):
        if P_[i] < P_ref:
            ind = i
            break

    # Integrate from P_ref to TOA
    out2 = integrate.solve_ivp(hydrostatic_equation, [P_ref, P_[-1]], np.array([0.0]), t_eval=P_[ind:], args=args, rtol=1e-6)
    assert out2.success
    # Integrate from P_ref to BOA
    out1 = integrate.solve_ivp(hydrostatic_equation, [P_ref, P_[0]], np.array([0.0]), t_eval=P_[:ind][::-1], args=args, rtol=1e-6)
    assert out1.success

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
        Pressure in dynes/cm^2
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
        raise ValueError('"P" must by an np.ndarray')
    if not isinstance(T, np.ndarray):
        raise ValueError('"P" must by an np.ndarray')
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
