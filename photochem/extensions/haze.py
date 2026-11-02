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
def _number_flux_loss(number_flux, mass_flux, state, particle_density, charging_radius):
    """Return coagulation loss per downward cm, in particles/(cm^3 s)."""
    if number_flux <= 0.0 or mass_flux <= 0.0:
        return 0.0
    mass = mass_flux/number_flux
    radius = np.cbrt(3.0*mass/(4.0*np.pi*particle_density))
    if not np.isfinite(radius) or radius <= 0.0:
        return np.inf
    temperature, density, number_density, gravity, production = state
    speed = fall_velocity(gravity, radius, particle_density, density,
                          number_density, temperature)
    coefficient = coagulation_coefficient(radius, particle_density,
                                          number_density, temperature, charging_radius)
    number = number_flux/speed
    return coefficient*number*number


@nb.njit
def _implicit_number_flux(available, mass_flux, state, step, particle_density,
                          seed_mass, charging_radius, maximum_mass):
    """Solve F + step*loss(F) = available inside physical flux bounds.

    A safeguarded Newton iteration uses a finite-difference slope and a
    bracket. The radius cap supplies a lower flux bound, limiting coagulation
    at the cap; new primary particles can subsequently release that bound.
    """
    lower = mass_flux/maximum_mass
    upper = min(available, mass_flux/seed_mass)
    if not np.isfinite(available):
        return 0.0, False
    if available < lower:
        # Project constrained stages onto the cap. A source tapering to zero
        # can otherwise drive an intermediate stage below the bound at every
        # step size, preventing the integrator from reaching an interval edge.
        return (lower, True) if lower > 0.0 else (0.0, False)
    if upper <= 0.0:
        return 0.0, True
    floor_loss = _number_flux_loss(lower, mass_flux, state, particle_density, charging_radius)
    if lower + step*floor_loss >= available:
        return lower, True
    # The sink is monotone in F for the adopted settling and coagulation laws.
    value = upper
    for iteration in range(100):
        loss = _number_flux_loss(value, mass_flux, state, particle_density, charging_radius)
        residual = value + step*loss - available
        if abs(residual) <= 2e-13*max(available, np.finfo(np.float64).tiny):
            return value, True
        if residual > 0.0 or not np.isfinite(residual):
            upper = value
        else:
            lower = value
        if upper - lower <= 2e-13*max(upper, np.finfo(np.float64).tiny):
            return .5*(lower + upper), True
        trial = .5*(lower + upper)
        if np.isfinite(loss) and value > 0.0:
            nearby = value*(1.0 - 1e-5)
            other_loss = _number_flux_loss(nearby, mass_flux, state, particle_density, charging_radius)
            slope = 1.0 + step*(loss - other_loss)/(value - nearby)
            newton = value - residual/slope
            if lower < newton < upper:
                trial = newton
        value = trial
    return value, False


@nb.njit
def _interval_mass_flux(distance, thickness, above, below, flux_above):
    """Analytically integrate the linearly interpolated mass source."""
    fraction = distance/thickness
    return flux_above + distance*(above[4] + .5*fraction*(below[4] - above[4]))


@nb.njit
def _number_flux_step(distance, number_flux, step, thickness, above, below,
                      flux_above, particle_density, seed_mass, charging_radius,
                      maximum_mass):
    """One second-order L-stable SDIRK step with exact source integration.

    Applying the two-stage method to the cumulative coagulation loss leaves
    mass production analytic. Negative unconstrained stages reject the step;
    capped stages are projected onto their minimum number flux. The diagonal
    coefficient is 1 - 1/sqrt(2).
    """
    gamma = 1.0 - 1.0/np.sqrt(2.0)
    middle = distance + gamma*step
    endpoint = distance + step
    c0 = _interval_mass_flux(distance, thickness, above, below, flux_above)
    c1 = _interval_mass_flux(middle, thickness, above, below, flux_above)
    c2 = _interval_mass_flux(endpoint, thickness, above, below, flux_above)
    state1 = (1.0 - middle/thickness)*above + middle/thickness*below
    state2 = (1.0 - endpoint/thickness)*above + endpoint/thickness*below
    available1 = number_flux + (c1 - c0)/seed_mass
    first, ok = _implicit_number_flux(available1, c1, state1, gamma*step,
                                      particle_density, seed_mass, charging_radius,
                                      maximum_mass)
    if not ok:
        return number_flux, False
    # Effective loss includes any suppression needed to respect the radius cap.
    used = available1 - first
    available2 = number_flux + (c2 - c0)/seed_mass - (1.0 - gamma)/gamma*used
    return _implicit_number_flux(available2, c2, state2, gamma*step,
                                  particle_density, seed_mass, charging_radius,
                                  maximum_mass)


@nb.njit
def _integrate_number_flux(number_flux, thickness, state_above, state_below,
                           flux_above, particle_density, seed_mass, charging_radius,
                           rtol, atol, maximum_mass=np.inf, max_steps=100000):
    """Advance downward number flux using adaptive implicit step doubling.

    Flux is in particles/(cm^2 s), thickness in cm, and endpoint states are
    [temperature K, gas density g/cm^3, gas number density cm^-3, gravity
    cm/s^2, mass production g/(cm^3 s)]. Mass flux is in g/(cm^2 s).
    Relative/absolute tolerances apply to log(F*m0/C), the surviving fraction
    of the primary-particle number flux. This controls fractional radius and
    number-flux error even when number production spans many orders of magnitude.
    Returns (number_flux, status): 0 success, 1 invalid/nonfinite state,
    2 minimum step, 3 attempted-step budget exhausted.
    """
    if (not np.isfinite(number_flux) or number_flux < 0.0 or
            not np.isfinite(flux_above) or flux_above < 0.0):
        return number_flux, 1
    distance = 0.0
    step = thickness
    rtol = max(rtol, 100*np.finfo(np.float64).eps)
    for attempt in range(max_steps):
        if step < 10*(np.nextafter(distance, np.inf) - distance) or distance + step == distance:
            return number_flux, 2
        endpoint = min(distance + step, thickness)
        step = endpoint - distance
        full, ok1 = _number_flux_step(distance, number_flux, step, thickness,
                                      state_above, state_below, flux_above,
                                      particle_density, seed_mass, charging_radius, maximum_mass)
        half, ok2 = _number_flux_step(distance, number_flux, .5*step, thickness,
                                      state_above, state_below, flux_above,
                                      particle_density, seed_mass, charging_radius, maximum_mass)
        fine, ok3 = _number_flux_step(distance + .5*step, half, .5*step, thickness,
                                      state_above, state_below, flux_above,
                                      particle_density, seed_mass, charging_radius, maximum_mass)
        error = np.inf
        c = _interval_mass_flux(endpoint, thickness, state_above, state_below, flux_above)
        if ok1 and ok2 and ok3 and np.isfinite(full) and np.isfinite(fine):
            if c == 0.0:
                error = 0.0
            elif full > 0.0 and fine > 0.0:
                log_full = np.log(full) + np.log(seed_mass) - np.log(c)
                log_fine = np.log(fine) + np.log(seed_mass) - np.log(c)
                scale = atol + rtol*max(abs(log_full), abs(log_fine))
                error = abs(log_fine - log_full)/(3.0*scale)
        if error <= 1.0:
            number_flux = fine
            distance = endpoint
            if distance >= thickness:
                return number_flux, 0
            step *= 5.0 if error == 0.0 else min(5.0, max(.2, .9*error**(-1.0/3.0)))
        else:
            step *= max(.2, .9*error**(-1.0/3.0)) if np.isfinite(error) else .2
    return number_flux, 3


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
    grid, with production of new primary particles. This extends the
    coagulation/sedimentation model of McKay et al. (1989),
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
        Radius of newly produced primary particles in cm; default 1e-7.
        Their mass m0 sets the number source J = production_rate/m0.
    charging_radius : float, optional
        Positive charge-suppression radius in cm; default 9e-6.
        Use np.inf to disable charge suppression.
    column_production : float, optional
        Internal column production in g/(cm^2 s). If supplied, scale the
        production profile to this integral over altitude before sweeping.
        Must be finite and nonnegative;
        a positive value requires a profile with positive integrated production.
    max_radius : float, optional
        Maximum particle radius in cm. Must be finite and at least seed_radius.
        Coagulation is limited at the imposed cap, conserving mass. Fresh
        primary-particle production can reduce radius below the cap again.
    sweep_below : float, optional
        Finite altitude in cm below which haze number density and mass flux
        are zero. Integration stops at the lowest retained sample, and a point
        exactly at the threshold is retained.
        Radius below it is held at the lowest retained sample's radius, or at
        seed_radius if no samples remain. Represents prescribed haze removal.
    rtol : float, optional
        Dimensionless relative tolerance on log(F*m0/C), the surviving
        primary-particle number-flux fraction; default 1e-6.
    atol : float, optional
        Dimensionless absolute tolerance on log(F*m0/C); default 1e-10.

    Returns
    -------
    particle_radius : ndarray, shape (N,)
        Particle radius in cm, ordered from bottom to top. An empty atmosphere
        retains seed_radius as a placeholder, except where sweep_below holds
        radius constant below the lowest retained sample.
    particle_number_density : ndarray, shape (N,)
        Particle number density in particles/cm^3, ordered from bottom to top.
    haze_mass_flux : ndarray, shape (N,)
        Positive downward mass flux in g/(cm^2 s), ordered from bottom to top.
        Zero at the upper boundary and below sweep_below. Without sweeping,
        the bottom value is total internal column production.

    Notes
    -----
    Production creates primary particles with number source J=q/m0, extending
    McKay's mass-only production assumption. Incoming mass and number fluxes
    are zero. Atmospheric profiles and production are interpolated linearly;
    cumulative mass production is integrated analytically. A Numba-compiled
    adaptive implicit number-flux integration handles stiff coagulation.
    Representative radius follows from C=mF wherever haze is present.
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

    seed_mass = (4.0/3.0)*np.pi*particle_density*seed_radius**3
    if not np.isfinite(seed_mass) or seed_mass <= 0.0:
        raise ValueError('Primary-particle mass must be finite and positive.')
    if not np.isfinite(haze_mass_flux[0]/seed_mass):
        raise ValueError('Integrated primary-particle production must be finite.')
    maximum_mass = np.inf
    if np.isfinite(radius_cap):
        maximum_mass = (4.0/3.0)*np.pi*particle_density*radius_cap**3
        if not np.isfinite(maximum_mass):
            raise ValueError('Maximum particle mass must be finite.')

    particle_radius = np.full(altitude.size, seed_radius, dtype=np.float64)
    particle_number_density = np.zeros(altitude.size, dtype=np.float64)
    first_retained = np.searchsorted(altitude, sweep_altitude, side='left')
    if haze_mass_flux[0] == 0.0 or first_retained >= altitude.size:
        haze_mass_flux.fill(0.0)
        return particle_radius, particle_number_density, haze_mass_flux

    number_flux = 0.0
    for i in range(altitude.size - 2, first_retained - 1, -1):
        if haze_mass_flux[i] == 0.0:
            continue
        if radius_cap == seed_radius:
            number_flux = haze_mass_flux[i]/seed_mass
        else:
            dz = altitude[i + 1] - altitude[i]
            state_above = np.array([
                temperature[i + 1], air_density[i + 1], air_number_density[i + 1],
                gravity[i + 1], production_rate[i + 1]
            ])
            state_below = np.array([
                temperature[i], air_density[i], air_number_density[i],
                gravity[i], production_rate[i]
            ])
            number_flux, status = _integrate_number_flux(
                number_flux, dz, state_above, state_below, haze_mass_flux[i + 1],
                particle_density, seed_mass, charging_radius, rtol, atol, maximum_mass
            )
            if status != 0:
                raise RuntimeError(f'Number-flux integration failed in altitude interval {i} (status {status}).')
        if not np.isfinite(number_flux) or number_flux <= 0.0:
            raise RuntimeError(f'Invalid particle number flux at altitude index {i}.')
        radius = np.cbrt(3.0*(haze_mass_flux[i]/number_flux)/(4.0*np.pi*particle_density))
        particle_radius[i] = max(seed_radius, min(radius, radius_cap))
        speed = fall_velocity(gravity[i], particle_radius[i], particle_density,
                              air_density[i], air_number_density[i], temperature[i])
        particle_number_density[i] = number_flux/speed
        if not np.isfinite(particle_number_density[i]) or particle_number_density[i] <= 0.0:
            raise RuntimeError(f'Invalid particle number density at altitude index {i}.')

    if first_retained > 0:
        particle_radius[:first_retained] = particle_radius[first_retained]
        haze_mass_flux[:first_retained] = 0.0
    return particle_radius, particle_number_density, haze_mass_flux
