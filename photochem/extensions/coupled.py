import numpy as np

from tempfile import NamedTemporaryFile

from photochem.utils import climate
from photochem.utils._format import yaml
from photochem.clima import AdiabatClimate, ClimaException
from photochem import EvoAtmosphere, PhotoException, PressureBC

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

        self.pc.var.nsteps_before_giveup = 7000
        self.pc.var.conv_longdy = 1e-3
        self.chemistry_converged = None
        self.chemistry_diagnostics = None
        self._edd = None
        self._edd_pressure = None
        self._climate_surface_partial_pressures = None
        self._solver_initialized = False
        self.niterations = 0
        self.converged = False
        self.give_up = False
        self.fixed_point_x = None
        self.fixed_point_gx = None
        self.temperature_residual = None
        self.log10pressure_residual = None
        self.solver_history = []

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

    def _surface_partial_pressures_from_photochem(self):
        """Return current fixed-pressure gas lower boundaries in dyn/cm².

        Other boundary types and particles are omitted. This reads the live
        model settings, including changes made through ``pc.set_lower_bc``.
        No initialized atmosphere is required.
        """
        pc = self.pc
        boundary_types = pc.var.lowerboundcond
        pressures = pc.var.lower_fix_press
        species_names = pc.dat.species_names
        return {
            species_names[i]: float(pressures[i])
            for i in range(pc.dat.np, pc.dat.nq)
            if boundary_types[i] == PressureBC
        }

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

    def _photochemical_profiles_from_climate(self, edd, edd_pressure, target_pressure):
        """Build independent photochemical input profiles from the climate.

        The caller validates Kzz and the target pressure before solving climate.
        Keep climate samples and Kzz knots within the photochemical domain,
        extending temperature and composition constantly above climate layers.
        Interpolate Kzz in log-pressure/log-Kzz space with constant endpoints.
        Neither model is modified.

        Parameters
        ----------
        edd : numpy.ndarray
            Validated scalar (zero-dimensional array) or Kzz profile in cm²/s.
        edd_pressure : numpy.ndarray or None
            Validated pressure knots in dyn/cm² for a Kzz profile.
        target_pressure : float
            Finite, positive photochemical top pressure in dyn/cm².

        Returns
        -------
        pressure, temperature, edd, mix : tuple
            Descending pressures (dyn/cm²), temperatures (K), Kzz (cm²/s),
            and a dictionary of gas volume-mixing-ratio profiles, ready for
            ``pc.initialize_atmosphere_p``.
        """
        pressure, temperature, mix = self.profiles_from_climate()
        if target_pressure >= pressure[0]:
            raise ValueError('Photochemical target_pressure must be below the climate surface pressure.')

        # Preserve climate samples and all Kzz knots inside the photochemical
        # domain. Include the climate top and the exact photochemical endpoint.
        grid_parts = [pressure[pressure > target_pressure], [target_pressure]]
        if target_pressure < self.c.P_top < pressure[0]:
            grid_parts.append([self.c.P_top])
        if edd.ndim == 1:
            inside = (edd_pressure <= pressure[0]) & (edd_pressure >= target_pressure)
            grid_parts.append(edd_pressure[inside])
        pressure_grid = np.unique(np.concatenate(grid_parts))[::-1].copy()
        log_pressure_grid = np.log10(pressure_grid)
        log_pressure_climate = np.log10(pressure[::-1])
        temperature_grid = np.interp(
            log_pressure_grid, log_pressure_climate, temperature[::-1]
        )
        mix_grid = {
            sp: np.interp(log_pressure_grid, log_pressure_climate, values[::-1])
            for sp, values in mix.items()
        }
        if edd.ndim == 0:
            edd_grid = np.full(pressure_grid.size, float(edd))
        else:
            edd_grid = 10.0**np.interp(
                log_pressure_grid, np.log10(edd_pressure[::-1]), np.log10(edd[::-1])
            )

        return pressure_grid, temperature_grid, edd_grid, mix_grid

    def initial_guess(self, edd, edd_pressure=None, surface_partial_pressures=None,
                      target_pressure=0.1, **kwargs):
        """Solve climate and initialize photochemistry without integrating it.

        Above the climate layers, temperature and gas mixing ratios are held
        constant. The photochemical grid extends to the TOA maintenance target.

        Parameters
        ----------
        edd : float or array-like
            Positive eddy diffusion in cm²/s. A scalar gives a constant
            profile; a vector requires ``edd_pressure``.
        edd_pressure : array-like, optional
            Finite, positive, strictly decreasing pressures in dyn/cm²,
            matching ``edd``. Kzz is interpolated in log-pressure/log-Kzz
            space and held constant outside the supplied range.
        surface_partial_pressures : dict[str, float], optional
            Initial climate surface partial pressures in dyn/cm². Defaults
            to the current fixed-pressure photochemical lower boundaries.
            Unspecified climate gases receive trace pressures of 1e-15.
        target_pressure : float, optional
            Photochemical top pressure in dyn/cm², by default 0.1. Sets both
            the initialization endpoint and TOA maintenance target.
        **kwargs
            Additional ``pc.initialize_atmosphere_p`` options, such as
            ``mode`` and ``tropopause_pressure``.

        Returns
        -------
        numpy.ndarray
            Climate fixed-point vector from ``fixed_point_state_from_climate``.
            The initialized photochemical atmosphere is not chemically converged.
        """

        # Validate eddy diffusion
        edd = np.asarray(edd, dtype=float)
        if edd.ndim > 1 or not np.all(np.isfinite(edd)) or np.any(edd <= 0.0):
            raise ValueError('edd must be a finite, positive scalar or one-dimensional profile.')
        if edd.ndim == 0:
            if edd_pressure is not None:
                raise ValueError('edd_pressure must be omitted when edd is a scalar.')
        else:
            if edd_pressure is None:
                raise ValueError('edd_pressure is required when edd is a profile.')
            edd_pressure = np.asarray(edd_pressure, dtype=float)
            if edd.size < 2 or edd_pressure.shape != edd.shape:
                raise ValueError('edd and edd_pressure must be matching one-dimensional arrays with at least two points.')
            if not np.all(np.isfinite(edd_pressure)) or np.any(edd_pressure <= 0.0):
                raise ValueError('edd_pressure must contain finite, positive pressures (dyn/cm²).')
            if np.any(np.diff(edd_pressure) >= 0.0):
                raise ValueError('edd_pressure must be strictly decreasing.')

        # Validate the target TOA pressure
        target_pressure = float(target_pressure)
        if not np.isfinite(target_pressure) or target_pressure <= 0.0:
            raise ValueError('Photochemical target_pressure must be finite and positive (dyn/cm²).')

        # Use the surface pressures to construct a sensible initial atmosphere
        if surface_partial_pressures is None:
            surface_partial_pressures = self._surface_partial_pressures_from_photochem()
        unknown_species = set(surface_partial_pressures) - set(self.c.species_names)
        if unknown_species:
            raise ValueError('Unknown climate species in surface_partial_pressures: '
                             + ', '.join(sorted(unknown_species)))
        if not surface_partial_pressures:
            raise ValueError(
                'No surface partial pressures are available for initialization. '
                'If the photochemical model has no fixed-pressure gas lower boundaries, '
                'you must supply surface_partial_pressures to initial_guess.'
            )
        P_i = np.full(len(self.c.species_names), 1e-15)
        for i, sp in enumerate(self.c.species_names):
            if sp in surface_partial_pressures:
                value = float(surface_partial_pressures[sp])
                if not np.isfinite(value) or value < 0.0:
                    raise ValueError(f'Surface partial pressure for {sp} must be finite and nonnegative (dyn/cm²).')
                P_i[i] = value
        if not any(float(value) > 0.0 for value in surface_partial_pressures.values()):
            raise ValueError('At least one supplied surface partial pressure must be positive.')

        # Reinitializing the atmosphere invalidates any previous solver session.
        self._solver_initialized = False
        self.converged = False
        self.give_up = False

        # Solve for steady-state climate.
        if not self.c.RCE_robust(P_i):
            raise ClimaException('Climate did not converge while constructing the initial guess.')

        # Extract the fixed point vector
        x = self.fixed_point_state_from_climate()

        # Build profiles covering the photochemical domain.
        pressure_grid, temperature_grid, edd_grid, mix_grid = (
            self._photochemical_profiles_from_climate(edd, edd_pressure, target_pressure)
        )

        self.pc.initialize_atmosphere_p(
            pressure=pressure_grid,
            temperature=temperature_grid,
            edd=edd_grid,
            mix=mix_grid,
            particle_radius=None,
            persistent=True,
            maintain_toa_pressure=True,
            target_pressure=target_pressure,
            **kwargs
        )
        # Retain independent prescribed inputs for subsequent map evaluations.
        self._edd = edd.copy()
        self._edd_pressure = None if edd_pressure is None else edd_pressure.copy()
        self._climate_surface_partial_pressures = {
            sp: float(P_i[i]) for i, sp in enumerate(self.c.species_names)
            if sp in self.climate_condensates
        }
        self.chemistry_converged = None
        self.chemistry_diagnostics = None
        return x

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

    def initialize_solver(self, edd, edd_pressure=None, surface_partial_pressures=None,
                          target_pressure=0.1, *, temperature_tol=1e-3,
                          log10pressure_tol=1e-3, relaxation=0.5, max_iterations=30,
                          require_chemistry_to_converge=False, **kwargs):
        """Initialize both the atmosphere and an incremental fixed-point solve.

        Calls ``initial_guess`` to solve the initial climate and initialize
        photochemistry. Chemical integration begins with the first ``iterate``.
        Restarting this solver resets its counters, residuals, and history.

        Parameters
        ----------
        edd : float or array-like
            Eddy diffusion input described by ``initial_guess``.
        edd_pressure : array-like, optional
            Pressure coordinates for a Kzz profile, as in ``initial_guess``.
        surface_partial_pressures : dict[str, float], optional
            Initial climate partial pressures, as in ``initial_guess``.
        target_pressure : float, optional
            Photochemical top pressure in dyn/cm², by default 0.1.
        temperature_tol : float, optional
            Maximum relative temperature mismatch, by default 1e-3.
        log10pressure_tol : float, optional
            Maximum absolute difference in log10 surface pressure (dex),
            by default 1e-3. Both tolerances must be finite and nonnegative.
        relaxation : float, optional
            Fraction of the unrelaxed update used for the next trial, in
            (0, 1], by default 0.5.
        max_iterations : int, optional
            Maximum number of completed map evaluations, by default 30.
        require_chemistry_to_converge : bool, optional
            Passed to ``fixed_point_map``, by default false. Coupled and
            chemical convergence are reported separately.
        **kwargs
            Additional options passed to ``initial_guess`` and its
            photochemical initializer, such as ``mode``.
        """
        temperature_tol = float(temperature_tol)
        log10pressure_tol = float(log10pressure_tol)
        relaxation = float(relaxation)
        if not np.isfinite(temperature_tol) or temperature_tol < 0.0:
            raise ValueError('temperature_tol must be finite and nonnegative.')
        if not np.isfinite(log10pressure_tol) or log10pressure_tol < 0.0:
            raise ValueError('log10pressure_tol must be finite and nonnegative (dex).')
        if not np.isfinite(relaxation) or not 0.0 < relaxation <= 1.0:
            raise ValueError('relaxation must be finite and in (0, 1].')
        if (isinstance(max_iterations, (bool, np.bool_)) or
                not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1):
            raise ValueError('max_iterations must be a positive integer.')

        self._solver_initialized = False
        x = self.initial_guess(
            edd, edd_pressure=edd_pressure,
            surface_partial_pressures=surface_partial_pressures,
            target_pressure=target_pressure, **kwargs
        )
        self.temperature_tol = temperature_tol
        self.log10pressure_tol = log10pressure_tol
        self.relaxation = relaxation
        self.max_iterations = int(max_iterations)
        self.require_chemistry_to_converge = bool(require_chemistry_to_converge)
        self._solver_next_x = x.copy()
        self.fixed_point_x = x.copy()
        self.fixed_point_gx = None
        self.niterations = 0
        self.converged = False
        self.give_up = False
        self.temperature_residual = None
        self.log10pressure_residual = None
        self.solver_history = []
        self._solver_initialized = True

    def iterate(self):
        """Perform one chemistry-to-climate map evaluation.

        Call ``initialize_solver`` first. Convergence uses the unrelaxed
        temperature and log-pressure mismatch. ``fixed_point_x`` stores the
        evaluated input, ``fixed_point_gx`` its climate result, and residuals
        and chemical status are recorded in ``solver_history``. The relaxed
        next trial is retained internally and is not applied to either model.

        Once converged or out of iterations, repeated calls return the same
        flags without further model evaluations. Model errors propagate and
        invalidate the solver session; reinitialize before continuing.

        Returns
        -------
        give_up : bool
            The iteration limit was reached without coupled convergence.
        converged : bool
            Both fixed-point residuals satisfy tolerance. This does not imply
            chemical convergence unless that was explicitly required; inspect
            ``chemistry_converged`` separately. ``solution`` remains the
            photochemical atmosphere from the last map evaluation.
        """
        if not self._solver_initialized:
            raise RuntimeError('Call initialize_solver before iterate, or reinitialize after a failed evaluation.')
        if self.converged or self.give_up:
            return self.give_up, self.converged

        x = self._solver_next_x.copy()
        # A failed model evaluation may have changed physical state. Do not
        # permit a retry using stale solver bookkeeping after an exception.
        self._solver_initialized = False
        gx = np.asarray(self.fixed_point_map(
            x, require_chemistry_to_converge=self.require_chemistry_to_converge
        ), dtype=float)
        if gx.shape != x.shape or not np.all(np.isfinite(gx)) or np.any(gx[1:] <= 0.0):
            raise ValueError('fixed_point_map returned an invalid pressure-temperature vector.')

        self.fixed_point_x = x.copy()
        self.fixed_point_gx = gx.copy()
        self.niterations += 1
        self.temperature_residual = float(np.max(np.abs(gx[1:] - x[1:])/x[1:]))
        self.log10pressure_residual = float(abs(gx[0] - x[0]))
        self.converged = (self.temperature_residual <= self.temperature_tol and
                          self.log10pressure_residual <= self.log10pressure_tol)
        self.give_up = not self.converged and self.niterations >= self.max_iterations
        self.solver_history.append({
            'iteration': self.niterations,
            'temperature_residual': self.temperature_residual,
            'log10pressure_residual': self.log10pressure_residual,
            'chemistry_converged': self.chemistry_converged,
        })
        if not self.converged and not self.give_up:
            self._solver_next_x = x + self.relaxation*(gx - x)
        self._solver_initialized = True
        return self.give_up, self.converged

    def find_steady_state(self, edd, **kwargs):
        """Initialize and iterate until coupled convergence or the iteration limit.

        Parameters
        ----------
        edd : float or array-like
            Eddy diffusion input for ``initial_guess``.
        **kwargs
            Atmospheric inputs and solver controls passed to
            ``initialize_solver``.

        Returns
        -------
        bool
            Whether the fixed-point solve converged. Inspect
            ``chemistry_converged`` separately. Counters, residuals, and the
            last photochemical atmosphere remain available after return.
        """
        self.initialize_solver(edd=edd, **kwargs)
        while True:
            give_up, converged = self.iterate()
            if give_up or converged:
                return converged

    def _climate_inputs_from_photochem(self):
        """Build climate pressures and dry-gas profiles from photochemistry.

        Estimate surface gas fractions by extending the lowest model layer
        constantly to the surface. Use the column-derived surface pressure
        and retain condensable inventories supplied to ``initial_guess``.
        Negative numerical traces are floored for climate interpolation only;
        the photochemical state is unchanged. Dry profiles are normalized
        within the dry gases, as required by ``AdiabatClimate.RCE``.
        """
        sol = self.solution()
        pressure = sol['pressure']
        if (not np.all(np.isfinite(pressure)) or np.any(pressure <= 0.0) or
                np.any(np.diff(pressure) >= 0.0)):
            raise PhotoException('Climate exchange requires finite, positive, strictly decreasing gas pressures.')
        if not np.all(np.isfinite(sol['density'])) or np.any(sol['density'] <= 0.0):
            raise PhotoException('Climate exchange requires finite, positive gas densities.')
        surface_pressure = self.pc.wrk.surface_pressure*1e6  # bars to dyn/cm²
        if not np.isfinite(surface_pressure) or surface_pressure <= 0.0:
            raise PhotoException('Photochemical surface pressure must be finite and positive.')

        gas_mix = np.array([sol[sp] for sp in self.c.species_names])
        if not np.all(np.isfinite(gas_mix)) or np.any(np.sum(gas_mix, axis=0) <= 0.0):
            raise PhotoException('Photochemical gas mixing ratios must be finite with a positive total.')
        gas_mix = np.maximum(gas_mix, 1e-30)
        gas_mix /= np.sum(gas_mix, axis=0)
        dry = np.array([sp not in self.climate_condensates for sp in self.c.species_names])
        P_i = np.empty(len(self.c.species_names))
        for i, sp in enumerate(self.c.species_names):
            if dry[i]:
                P_i[i] = gas_mix[i, 0]*surface_pressure
            else:
                P_i[i] = self._climate_surface_partial_pressures[sp]

        custom_dry_mix = None
        if np.any(dry):
            dry_mix = gas_mix[dry]/np.sum(gas_mix[dry], axis=0)
            custom_dry_mix = {'pressure': pressure.copy()}
            for sp, values in zip(np.asarray(self.c.species_names)[dry], dry_mix):
                custom_dry_mix[str(sp)] = values.copy()
        return P_i, custom_dry_mix

    def fixed_point_map(self, x, require_chemistry_to_converge=False):
        """Evaluate the chemistry-to-climate map without advancing a solver.

        The input surface pressure defines temperature coordinates; it is not
        imposed on the photochemical atmosphere. Kzz and condensable inventories
        come from the last successful ``initial_guess``. Chemistry uses its
        existing composition as a warm start, while climate must converge.
        The returned climate state is not applied back to photochemistry.

        Parameters
        ----------
        x : array-like
            ``[log10(P_surface), T_surface, ..., T_top]`` in the representation
            returned by ``fixed_point_state_from_climate``.
        require_chemistry_to_converge : bool, optional
            If false (default), chemistry may exhaust its accepted-step budget
            and still supply climate inputs. Other integration failures remain
            errors. If true, chemical nonconvergence always raises an error.

        Returns
        -------
        numpy.ndarray
            Climate-derived fixed-point vector. ``chemistry_converged`` and
            ``chemistry_diagnostics`` describe this evaluation's chemical run;
            neither indicates convergence of a coupled fixed-point solver.
        """
        if self._edd is None or not self.pc.atmosphere_initialized:
            raise PhotoException('Call initial_guess before evaluating fixed_point_map.')
        x = np.asarray(x, dtype=float)
        if x.shape != (self.c.T.size + 2,) or not np.all(np.isfinite(x)):
            raise ValueError('x must be a finite vector containing log10 surface pressure and one temperature per coupling-grid point.')
        if np.any(x[1:] <= 0.0):
            raise ValueError('Fixed-point temperatures must be positive (K).')
        log_pressure = self.make_pressure_grid(x[0])
        if x[0] > np.log10(np.finfo(float).max):
            raise ValueError('Fixed-point surface pressure exceeds the floating-point range.')
        pressure = 10.0**log_pressure
        pc, c = self.pc, self.c
        profile = pc.var.press_temp_edd_profile
        maintenance = pc.var.toa_pressure_maintenance
        if not profile.enabled:
            raise PhotoException('fixed_point_map requires an enabled pressure-temperature-eddy profile.')
        target_pressure = (maintenance.target_pressure if maintenance.enabled
                           else profile.pressure[-1])

        # Include the actual surface and original Kzz knots without imposing
        # x's surface pressure or rebuilding the chemical composition.
        bottom_pressure = max(pressure[0], pc.wrk.surface_pressure*1e6)
        if (not np.isfinite(bottom_pressure) or not np.isfinite(target_pressure) or
                target_pressure <= 0.0 or target_pressure >= bottom_pressure):
            raise ValueError('Photochemical pressure endpoints must be finite, positive, and strictly decreasing.')
        parts = [pressure[pressure > target_pressure], [bottom_pressure, target_pressure]]
        if self._edd_pressure is not None:
            inside = ((self._edd_pressure <= bottom_pressure) &
                      (self._edd_pressure >= target_pressure))
            parts.append(self._edd_pressure[inside])
        pressure_grid = np.unique(np.concatenate(parts))[::-1].copy()
        temperature_grid = np.interp(
            np.log10(pressure_grid), log_pressure[::-1], x[1:][::-1]
        )
        if self._edd.ndim == 0:
            edd_grid = np.full(pressure_grid.size, float(self._edd))
        else:
            edd_grid = 10.0**np.interp(
                np.log10(pressure_grid), np.log10(self._edd_pressure[::-1]),
                np.log10(self._edd[::-1])
            )
        self.chemistry_converged = None
        self.chemistry_diagnostics = None
        pc.destroy_stepper()
        pc.set_press_temp_edd_profile(
            pressure_grid, temperature_grid, edd_grid,
            trop_p=profile.trop_p, mode=profile.mode,
            hydro_pressure=profile.hydro_pressure,
            maintain_toa_pressure=maintenance.enabled, target_pressure=target_pressure
        )
        try:
            converged = pc.find_steady_state()
            self.chemistry_converged = bool(converged)
            self.chemistry_diagnostics = {
                'converged': bool(converged),
                'steps': pc.wrk.nsteps_total,
                'integration_errors': pc.wrk.nerrors_total,
                'resync_failures': pc.wrk.n_resync_failures,
                'longdy': pc.wrk.longdy,
                'longdydt': pc.wrk.longdydt,
                'time': pc.wrk.tn,
                'step_limit_reached': pc.wrk.nsteps_total >= pc.var.nsteps_before_giveup,
            }
        finally:
            pc.destroy_stepper()
        if not converged:
            if require_chemistry_to_converge:
                raise PhotoException('Photochemistry did not converge during fixed_point_map.')
            if not self.chemistry_diagnostics['step_limit_reached']:
                raise PhotoException('Photochemistry gave up before reaching its step limit during fixed_point_map.')

        P_i, custom_dry_mix = self._climate_inputs_from_photochem()
        # Climate layers are at the midpoints of equally spaced log pressures.
        T_init = np.concatenate(([x[1]], 0.5*(x[1:-1] + x[2:])))
        converged = c.RCE_robust(
            P_i, custom_dry_mix, T_init=T_init,
            convecting_with_below_init=c.convecting_with_below
        )
        if not converged:
            converged = c.RCE_robust(P_i, custom_dry_mix)
        if not converged:
            raise ClimaException('Climate did not converge during fixed_point_map.')
        return self.fixed_point_state_from_climate()


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
