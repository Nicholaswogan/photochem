import numpy as np
import numba as nb
from scipy import constants

K_BOLTZ = constants.Boltzmann*1.0e7  # erg/K


@nb.njit
def gas_mean_free_path(air_number_density):
    """Estimate the gas mean free path using a fixed effective collision area.

    Parameters
    ----------
    air_number_density : float
        Gas number density in molecules/cm^3.

    Returns
    -------
    float
        Mean free path in cm, using an effective collision area of 6e-15 cm^2.
    """
    return 1.0/(air_number_density*6.0e-15)


@nb.njit
def slip_correction_factor(particle_radius, air_number_density):
    """Compute the Cunningham slip correction for a spherical particle.

    Parameters
    ----------
    particle_radius : float
        Particle radius in cm.
    air_number_density : float
        Gas number density in molecules/cm^3.

    Returns
    -------
    float
        Dimensionless factor multiplying the Stokes settling speed.

    References
    ----------
    Seinfeld and Pandis (2006), Atmospheric Chemistry and Physics, Eq. 9.34.
    """
    free_path = gas_mean_free_path(air_number_density)
    correction_factor = 1.0 + (free_path/particle_radius)*(
        1.257 + 0.4*np.exp(-1.1*particle_radius/free_path)
    )
    return correction_factor


@nb.njit
def dynamic_viscosity_air(temperature):
    """Compute Earth-air dynamic viscosity using Sutherland's relation.

    Parameters
    ----------
    temperature : float
        Gas temperature in K.

    Returns
    -------
    float
        Dynamic viscosity in g/(cm s).

    References
    ----------
    White (2006), Viscous Fluid Flow, Eq. 1-36 and Table 1-2 (air).
    """
    temperature_ref = 273.0
    viscosity_ref = 1.716e-4
    sutherland_temperature = 111.0
    return viscosity_ref*(temperature/temperature_ref)**1.5*(
        temperature_ref + sutherland_temperature
    )/(temperature + sutherland_temperature)


@nb.njit
def stokes_law(gravity, particle_radius, particle_density, air_density, viscosity):
    """Compute the signed downward Stokes terminal speed of a sphere.

    Parameters
    ----------
    gravity : float
        Magnitude of gravitational acceleration in cm/s^2.
    particle_radius : float
        Particle radius in cm.
    particle_density : float
        Particle material density in g/cm^3.
    air_density : float
        Gas mass density in g/cm^3.
    viscosity : float
        Gas dynamic viscosity in g/(cm s).

    Returns
    -------
    float
        Settling velocity in cm/s, positive downward for particles denser
        than the gas. Assumes low-Reynolds-number drag.

    References
    ----------
    Seinfeld and Pandis (2006), Atmospheric Chemistry and Physics, Eq. 9.29.
    """
    return (2.0/9.0)*gravity*particle_radius**2*(
        particle_density - air_density
    )/viscosity


@nb.njit
def fall_velocity(gravity, particle_radius, particle_density, air_density,
                  air_number_density, temperature):
    """Compute slip-corrected settling using Earth-air viscosity.

    Parameters
    ----------
    gravity : float
        Magnitude of gravitational acceleration in cm/s^2.
    particle_radius : float
        Particle radius in cm.
    particle_density : float
        Particle material density in g/cm^3.
    air_density : float
        Gas mass density in g/cm^3.
    air_number_density : float
        Gas number density in molecules/cm^3, used for the slip correction.
    temperature : float
        Gas temperature in K.

    Returns
    -------
    float
        Terminal settling velocity in cm/s, positive downward for particles
        denser than the gas.

    References
    ----------
    Seinfeld and Pandis (2006), Atmospheric Chemistry and Physics,
    Eqs. 9.29 and 9.34; White (2006), Viscous Fluid Flow, Eq. 1-36.
    """
    viscosity = dynamic_viscosity_air(temperature)
    stokes_velocity = stokes_law(
        gravity, particle_radius, particle_density, air_density, viscosity
    )
    correction_factor = slip_correction_factor(particle_radius, air_number_density)
    return stokes_velocity*correction_factor


@nb.njit
def coagulation_kernel_free_molecular(particle_radius, particle_density, temperature):
    """Compute the free-molecular Brownian encounter kernel for equal spheres.

    Parameters
    ----------
    particle_radius : float
        Particle radius in cm.
    particle_density : float
        Particle material density in g/cm^3.
    temperature : float
        Gas and particle temperature in K.

    Returns
    -------
    float
        Uncharged pair-encounter kernel in cm^3/s, before division by two
        to obtain the particle-number loss coefficient.
    """
    return 8.0*np.sqrt(3.0*K_BOLTZ*temperature/particle_density)*np.sqrt(particle_radius)


@nb.njit
def coagulation_kernel_continuum(temperature, viscosity):
    """Compute the continuum Brownian encounter kernel for equal spheres.

    Parameters
    ----------
    temperature : float
        Gas temperature in K.
    viscosity : float
        Gas dynamic viscosity in g/(cm s).

    Returns
    -------
    float
        Uncharged pair-encounter kernel in cm^3/s, before division by two
        to obtain the particle-number loss coefficient.
    """
    return 8.0*K_BOLTZ*temperature/(3.0*viscosity)


@nb.njit
def knudsen_number(air_number_density, particle_radius):
    """Compute the gas mean free path divided by the particle radius.

    Parameters
    ----------
    air_number_density : float
        Gas number density in molecules/cm^3.
    particle_radius : float
        Particle radius in cm.

    Returns
    -------
    float
        Dimensionless Knudsen number, using particle radius as the length scale.
    """
    return gas_mean_free_path(air_number_density)/particle_radius


@nb.njit
def coagulation_kernel(particle_radius, particle_density, air_number_density,
                       temperature):
    """Blend the free-molecular and continuum Brownian encounter kernels.

    Parameters
    ----------
    particle_radius : float
        Particle radius in cm.
    particle_density : float
        Particle material density in g/cm^3.
    air_number_density : float
        Gas number density in molecules/cm^3.
    temperature : float
        Gas and particle temperature in K; sets the Earth-air viscosity.

    Returns
    -------
    float
        Uncharged pair-encounter kernel in cm^3/s, before charge suppression
        or division by two to obtain the particle-number loss coefficient.
    """
    viscosity = dynamic_viscosity_air(temperature)
    kernel_free = coagulation_kernel_free_molecular(
        particle_radius, particle_density, temperature
    )
    kernel_continuum = coagulation_kernel_continuum(temperature, viscosity)
    kn = knudsen_number(air_number_density, particle_radius)
    return (kn*kernel_free + kernel_continuum)/(1.0 + kn)


@nb.njit
def coagulation_coefficient(particle_radius, particle_density, air_number_density,
                            temperature, charging_radius):
    """Compute the charge-suppressed particle-number loss coefficient.

    Parameters
    ----------
    particle_radius : float
        Particle radius in cm.
    particle_density : float
        Particle material density in g/cm^3.
    air_number_density : float
        Gas number density in molecules/cm^3.
    temperature : float
        Gas and particle temperature in K; sets the Earth-air viscosity.
    charging_radius : float
        Positive charge-suppression radius in cm; the rate is reduced by
        a factor of e when particle_radius equals charging_radius.

    Returns
    -------
    float
        Coefficient K in cm^3/s for the particle-number sink -K*n^2.
        Includes exp(-particle_radius/charging_radius) and division by two
        to avoid counting each particle pair twice.

    References
    ----------
    McKay, Pollack, and Courtin (1989), The thermal structure of Titan's
    atmosphere, Icarus 80, 23-53, doi:10.1016/0019-1035(89)90160-7.
    """
    kernel = coagulation_kernel(
        particle_radius, particle_density, air_number_density, temperature
    )
    return 0.5*kernel*np.exp(-particle_radius/charging_radius)


@nb.njit(error_model="numpy")
def _mckay_log_radius_rhs(distance, log_radius, thickness, state_above,
                          state_below, flux_above, particle_density,
                          charging_radius):
    """Return d(log r)/d(downward distance), in cm^-1, within one interval."""
    fraction = distance/thickness
    state = (1.0 - fraction)*state_above + fraction*state_below
    temperature, density, number_density, gravity, production = state
    flux = flux_above + distance*(
        state_above[4] + 0.5*fraction*(state_below[4] - state_above[4])
    )
    if not np.isfinite(log_radius):
        return np.nan
    radius = np.exp(log_radius)
    if not np.isfinite(radius) or radius <= 0.0:
        return np.nan
    speed = fall_velocity(
        gravity, radius, particle_density, density, number_density, temperature
    )
    coefficient = coagulation_coefficient(
        radius, particle_density, number_density, temperature, charging_radius
    )
    gradient = production/(3.0*flux) + coefficient*flux/(
        4.0*np.pi*particle_density*radius**3*speed**2
    )
    return gradient


# Dormand-Prince 5(4) tableau; the seventh stage is the next step's first RHS.
_RK45_C = np.array([0.0, 1/5, 3/10, 4/5, 8/9, 1.0])
_RK45_A = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0],
    [1/5, 0.0, 0.0, 0.0, 0.0],
    [3/40, 9/40, 0.0, 0.0, 0.0],
    [44/45, -56/15, 32/9, 0.0, 0.0],
    [19372/6561, -25360/2187, 64448/6561, -212/729, 0.0],
    [9017/3168, -355/33, 46732/5247, 49/176, -5103/18656],
])
_RK45_B = np.array([35/384, 0.0, 500/1113, 125/192, -2187/6784, 11/84])
_RK45_E = np.array([-71/57600, 0.0, 71/16695, -71/1920, 17253/339200, -22/525, 1/40])


@nb.njit
def _integrate_log_radius(log_radius, thickness, state_above, state_below,
                          flux_above, particle_density, charging_radius,
                          rtol, atol, log_max_radius=np.inf, max_steps=100000):
    """Advance scalar log radius downward using adaptive Dormand-Prince RK45.

    Parameters
    ----------
    log_radius : float
        Initial ln(radius / 1 cm).
    thickness : float
        Positive interval thickness in cm.
    state_above, state_below : ndarray, shape (5,)
        Endpoint temperature (K), gas mass density (g/cm^3), gas number density
        (molecules/cm^3), gravity (cm/s^2), and production (g/(cm^3 s)).
    flux_above : float
        Positive downward mass flux at the upper endpoint in g/(cm^2 s).
    particle_density : float
        Particle material density in g/cm^3.
    charging_radius : float
        Charge-suppression radius in cm; np.inf disables suppression.
    rtol, atol : float
        Positive dimensionless relative and absolute log-radius tolerances.
    log_max_radius : float, optional
        ln(maximum radius / 1 cm); default np.inf disables the cap.
    max_steps : int, optional
        Maximum attempted steps, including rejections; default 100000.

    Returns
    -------
    log_radius : float
        Final log radius, or the last accepted value on failure.
    capped : bool
        Whether an accepted state reached the maximum radius.
    status : int
        0: success; 1: nonfinite RHS; 2: minimum step; 3: step budget exhausted.

    References
    ----------
    Dormand and Prince (1980), A family of embedded Runge-Kutta formulae,
    Journal of Computational and Applied Mathematics 6, 19-26.
    """
    if log_radius >= log_max_radius:
        return log_max_radius, True, 0
    gradient = _mckay_log_radius_rhs(
        0.0, log_radius, thickness, state_above, state_below,
        flux_above, particle_density, charging_radius
    )
    if not np.isfinite(gradient):
        return log_radius, False, 1

    # Keep the first trial near the seed small; subsequent steps use error control.
    step = min(thickness, 0.05/gradient) if gradient > 0.0 else thickness
    rtol = max(rtol, 100.0*np.finfo(np.float64).eps)
    distance = 0.0
    slopes = np.empty(7, dtype=np.float64)
    rejected = False
    nonfinite_trial = False
    for attempt in range(max_steps):
        minimum_step = 10.0*(np.nextafter(distance, np.inf) - distance)
        if step < minimum_step or distance + step == distance:
            return log_radius, False, 1 if nonfinite_trial else 2
        endpoint = min(distance + step, thickness)
        step = endpoint - distance

        slopes[0] = gradient
        finite = True
        for stage in range(1, 6):
            increment = 0.0
            for j in range(stage):
                increment += _RK45_A[stage, j]*slopes[j]
            slopes[stage] = _mckay_log_radius_rhs(
                distance + _RK45_C[stage]*step, log_radius + step*increment,
                thickness, state_above, state_below, flux_above,
                particle_density, charging_radius
            )
            if not np.isfinite(slopes[stage]):
                finite = False
                break

        error_norm = np.inf
        candidate = log_radius
        if finite:
            increment = 0.0
            for j in range(6):
                increment += _RK45_B[j]*slopes[j]
            candidate = log_radius + step*increment
            slopes[6] = _mckay_log_radius_rhs(
                endpoint, candidate, thickness, state_above, state_below,
                flux_above, particle_density, charging_radius
            )
            finite = np.isfinite(candidate) and np.isfinite(slopes[6])
            if finite:
                error = 0.0
                for j in range(7):
                    error += _RK45_E[j]*slopes[j]
                scale = atol + rtol*max(abs(log_radius), abs(candidate))
                error_norm = abs(step*error)/scale

        if np.isfinite(error_norm) and error_norm < 1.0:
            distance = endpoint
            log_radius = candidate
            gradient = slopes[6]
            if log_radius >= log_max_radius:
                return log_max_radius, True, 0
            if distance >= thickness:
                return log_radius, False, 0
            factor = 10.0 if error_norm == 0.0 else min(10.0, 0.9*error_norm**(-0.2))
            if rejected:
                factor = min(1.0, factor)
            step *= factor
            rejected = False
            nonfinite_trial = False
        else:
            nonfinite_trial = not finite or not np.isfinite(error_norm)
            factor = max(0.2, 0.9*error_norm**(-0.2)) if np.isfinite(error_norm) else 0.2
            step *= factor
            rejected = True

    return log_radius, False, 3


@nb.njit
def _validate_haze_profile(profile, size):
    """Validate the dimensions, length, and finite values of one point profile."""
    if profile.ndim != 1 or profile.size != size:
        raise ValueError('All input profiles must be one-dimensional and match altitude.')
    if not np.all(np.isfinite(profile)):
        raise ValueError('All input profiles must contain finite values.')


@nb.njit
def mckay_microphysics(
    altitude,
    temperature,
    air_density,
    air_number_density,
    gravity,
    production_rate,
    particle_density=1.0,
    seed_radius=1.0e-7,
    charging_radius=9.0e-6,
    column_production=None,
    max_radius=None,
    sweep_below=None,
    rtol=1.0e-6,
    atol=1.0e-10,
):
    """Solve steady monodisperse haze microphysics on a supplied altitude
    grid. The model generally follows the one used in McKay et al. (1989), 
    DOI: 10.1016/0019-1035(89)90160-7.

    Parameters
    ----------
    altitude : ndarray, shape (N,)
        Sample altitudes in cm, strictly increasing from bottom to top.
        These are point values, not cell edges for cell-averaged inputs.
    temperature : ndarray, shape (N,)
        Gas and particle temperature in K.
    air_density : ndarray, shape (N,)
        Gas mass density in g/cm^3, less than particle_density.
    air_number_density : ndarray, shape (N,)
        Gas number density in molecules/cm^3.
    gravity : ndarray, shape (N,)
        Positive gravitational acceleration magnitude in cm/s^2.
    production_rate : ndarray, shape (N,)
        Nonnegative haze mass production at each sample point in g/(cm^3 s).
    particle_density : float, optional
        Constant particle material density in g/cm^3; default 1.
    seed_radius : float, optional
        Radius in cm at the highest point with positive cumulative production;
        default 1e-7. This prescribes the initial particle population.
    charging_radius : float, optional
        Positive charge-suppression radius in cm; default 9e-6.
        Use np.inf to disable charge suppression.
    column_production : float, optional
        Total column production in g/(cm^2 s). If supplied, scale the production
        profile to this integral over altitude before any sweeping is applied.
        Must be finite and nonnegative;
        a positive value requires a profile with positive integrated production.
    max_radius : float, optional
        Maximum particle radius in cm. Must be finite and at least seed_radius.
        Integration stops growing particles when this imposed cap is reached;
        number density is still calculated from the prescribed mass flux.
    sweep_below : float, optional
        Finite altitude in cm below which haze number density and mass flux
        are zero. Integration stops at the lowest retained sample, and a point
        exactly at the threshold is retained.
        Radius below it is held at the lowest retained sample's radius, or at
        seed_radius if no samples remain. Represents prescribed haze removal.
    rtol : float, optional
        Positive relative tolerance for integrating log radius; default 1e-6.
        Dimensionless and applied to ln(radius / 1 cm), not directly to radius.
    atol : float, optional
        Positive absolute tolerance for integrating log radius; default 1e-10.
        Dimensionless; small errors in log radius are fractional radius errors.

    Returns
    -------
    particle_radius : ndarray, shape (N,)
        Particle radius in cm, ordered from bottom to top. Empty regions
        retain seed_radius as a placeholder, except where sweep_below holds
        radius constant below the lowest retained sample.
    particle_number_density : ndarray, shape (N,)
        Particle number density in particles/cm^3, ordered from bottom to top.
    haze_mass_flux : ndarray, shape (N,)
        Positive downward mass flux in g/(cm^2 s), ordered from bottom to top.
        Zero at the upper boundary and below sweep_below. Without sweeping,
        the bottom value is total column production.

    Notes
    -----
    Input validation, the numerical driver, and scalar adaptive
    Dormand-Prince RK45 integration are Numba compiled. Atmospheric profiles
    and production are interpolated linearly,
    and cumulative production is integrated analytically within each interval.
    Integration restarts at input points and stops once an accepted step
    reaches max_radius; its precise crossing altitude is not calculated.

    """

    # Check input floats
    particle_density = float(particle_density)
    seed_radius = float(seed_radius)
    charging_radius = float(charging_radius)
    rtol = float(rtol)
    atol = float(atol)
    if not np.isfinite(rtol) or rtol <= 0.0 or not np.isfinite(atol) or atol <= 0.0:
        raise ValueError('rtol and atol must be finite and positive.')
    column_target = -1.0
    if column_production is not None:
        column_target = float(column_production)
        if not np.isfinite(column_target) or column_target < 0.0:
            raise ValueError('column_production must be finite and nonnegative.')

    # Convert and check input arrays
    altitude = np.asarray(altitude, dtype=np.float64)
    temperature = np.asarray(temperature, dtype=np.float64)
    air_density = np.asarray(air_density, dtype=np.float64)
    air_number_density = np.asarray(air_number_density, dtype=np.float64)
    gravity = np.asarray(gravity, dtype=np.float64)
    production_rate = np.asarray(production_rate, dtype=np.float64)
    if altitude.ndim != 1 or altitude.size < 2:
        raise ValueError('altitude must be a one-dimensional array with at least two points.')
    altitude = altitude.ravel()
    _validate_haze_profile(temperature, altitude.size)
    _validate_haze_profile(air_density, altitude.size)
    _validate_haze_profile(air_number_density, altitude.size)
    _validate_haze_profile(gravity, altitude.size)
    _validate_haze_profile(production_rate, altitude.size)
    # Keep downstream array types one-dimensional even for rejected inputs.
    temperature = temperature.ravel()
    air_density = air_density.ravel()
    air_number_density = air_number_density.ravel()
    gravity = gravity.ravel()
    production_rate = production_rate.ravel()
    if not np.all(np.isfinite(altitude)) or np.any(np.diff(altitude) <= 0.0):
        raise ValueError('altitude must be finite and strictly increasing from bottom to top.')
    if not np.isfinite(particle_density) or particle_density <= 0.0:
        raise ValueError('particle_density must be finite and positive.')
    if not np.isfinite(seed_radius) or seed_radius <= 0.0:
        raise ValueError('seed_radius must be finite and positive.')
    if np.isnan(charging_radius) or charging_radius <= 0.0:
        raise ValueError('charging_radius must be positive (or np.inf).')
    radius_cap = np.inf
    if max_radius is not None:
        radius_cap = float(max_radius)
        if not np.isfinite(radius_cap) or radius_cap < seed_radius:
            raise ValueError('max_radius must be finite and at least seed_radius.')
    sweep_altitude = -np.inf
    if sweep_below is not None:
        sweep_altitude = float(sweep_below)
        if not np.isfinite(sweep_altitude):
            raise ValueError('sweep_below must be finite.')
    if np.any(temperature <= 0.0) or np.any(air_number_density <= 0.0) or np.any(gravity <= 0.0):
        raise ValueError('Temperature, gas number density, and gravity must be positive.')
    if np.any(air_density < 0.0) or np.any(air_density >= particle_density):
        raise ValueError('Gas mass density must be nonnegative and less than particle_density.')
    if np.any(production_rate < 0.0):
        raise ValueError('production_rate must be nonnegative.')

    # Compute the production between each altitude interval (g/cm^2/s)
    layer_production = (
        0.5*(production_rate[:-1] + production_rate[1:])*np.diff(altitude)
    )
    # Compute the cumulative production rate. It is zero at the TOA, and then
    # the total production rate (g/cm^2/s) at the first index.
    haze_mass_flux = np.zeros(altitude.size, dtype=np.float64)
    haze_mass_flux[:-1] = np.cumsum(layer_production[::-1])[::-1]
    if not np.all(np.isfinite(haze_mass_flux)):
        raise ValueError('Integrated haze production must be finite.')

    if column_target >= 0.0:
        # If a specific column production is specified, then we adjust the production
        # rate profile accordingly.
        if column_target == 0.0:
            # If no production, then we take this shortcut.
            production_rate = np.zeros_like(production_rate)
            haze_mass_flux.fill(0.0)
        elif haze_mass_flux[0] == 0.0:
            # Guards against a division by zero.
            raise ValueError('Positive column_production requires a nonzero production profile.')
        else:
            # Scale the production rate by the column production. Compute a new haze_mass_flux.
            production_rate = production_rate*(column_target/haze_mass_flux[0])
            layer_production = (
                0.5*(production_rate[:-1] + production_rate[1:])*np.diff(altitude)
            )
            haze_mass_flux[:-1] = np.cumsum(layer_production[::-1])[::-1]
            if (not np.all(np.isfinite(production_rate)) or
                    not np.all(np.isfinite(haze_mass_flux)) or haze_mass_flux[0] <= 0.0):
                raise ValueError('Normalized haze production must be finite and have a positive integral.')

    # Prepare the outputs. We start with the seed_radius then we will overwrite.
    particle_radius = np.full(altitude.size, seed_radius, dtype=np.float64)
    particle_number_density = np.zeros(altitude.size, dtype=np.float64)
    first_retained = np.searchsorted(altitude, sweep_altitude, side='left')

    # Determine the highest point in the atmosphere, where there is
    # non-zero haze production.
    seed_index = -1
    for i in range(altitude.size - 1, -1, -1):
        if haze_mass_flux[i] > 0.0:
            seed_index = i
            break
    if seed_index == -1 or seed_index < first_retained:
        # No haze remains above the sweep threshold.
        haze_mass_flux.fill(0.0)
        return particle_radius, particle_number_density, haze_mass_flux

    # Initialize the number density at the prescribed seed radius.
    speed = fall_velocity(
        gravity[seed_index], seed_radius, particle_density, air_density[seed_index],
        air_number_density[seed_index], temperature[seed_index]
    )
    mass = (4.0/3.0)*np.pi*particle_density*seed_radius**3
    mass_flux_per_particle = mass*speed
    if not np.isfinite(mass_flux_per_particle) or mass_flux_per_particle <= 0.0:
        raise RuntimeError(f'Invalid particle settling flux at altitude index {seed_index}.')
    particle_number_density[seed_index] = haze_mass_flux[seed_index]/mass_flux_per_particle
    if (not np.isfinite(particle_number_density[seed_index]) or
            particle_number_density[seed_index] <= 0.0):
        raise RuntimeError(f'Invalid particle number density at altitude index {seed_index}.')

    log_max_radius = np.log(radius_cap)
    capped = radius_cap == seed_radius

    # From the seed point downward, we do a loop
    for i in range(seed_index - 1, first_retained - 1, -1):

        if capped:
            particle_radius[i] = radius_cap
        else:
            # Between each interval, integrate the radius equation through each interval.
            # In practice we integrate the natural log of radius, for numerical purposes.

            # Prepare inputs
            dz = altitude[i + 1] - altitude[i]
            state_above = np.array([
                temperature[i + 1], air_density[i + 1], air_number_density[i + 1],
                gravity[i + 1], production_rate[i + 1]
            ])
            state_below = np.array([
                temperature[i], air_density[i], air_number_density[i],
                gravity[i], production_rate[i]
            ])
            flux_above = haze_mass_flux[i + 1]
            log_radius = np.log(particle_radius[i + 1])

            # Run adaptive scalar RK45 through this interval.
            log_radius, capped, status = _integrate_log_radius(
                log_radius, dz, state_above, state_below, flux_above,
                particle_density, charging_radius, rtol, atol, log_max_radius
            )
            if status != 0:
                if status == 1:
                    raise RuntimeError(f'Nonfinite radius integration in altitude interval {i}.')
                if status == 2:
                    raise RuntimeError(f'Radius integration step too small in altitude interval {i}.')
                raise RuntimeError(f'Radius integration exceeded its step budget in altitude interval {i}.')
            particle_radius[i] = radius_cap if capped else np.exp(log_radius)

            if not np.isfinite(particle_radius[i]) or particle_radius[i] <= 0.0:
                raise RuntimeError(f'Invalid particle radius in altitude interval {i}.')

        # Compute the particle number density from the radius
        speed = fall_velocity(
            gravity[i], particle_radius[i], particle_density, air_density[i],
            air_number_density[i], temperature[i]
        )
        mass = (4.0/3.0)*np.pi*particle_density*particle_radius[i]**3
        mass_flux_per_particle = mass*speed
        if not np.isfinite(mass_flux_per_particle) or mass_flux_per_particle <= 0.0:
            raise RuntimeError(f'Invalid particle settling flux at altitude index {i}.')
        particle_number_density[i] = haze_mass_flux[i]/mass_flux_per_particle
        if not np.isfinite(particle_number_density[i]) or particle_number_density[i] <= 0.0:
            raise RuntimeError(f'Invalid particle number density at altitude index {i}.')

    if first_retained > 0:
        particle_radius[:first_retained] = particle_radius[first_retained]
        haze_mass_flux[:first_retained] = 0.0

    return particle_radius, particle_number_density, haze_mass_flux
