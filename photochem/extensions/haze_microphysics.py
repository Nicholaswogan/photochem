import numpy as np
import numba as nb
from scipy import constants
from scipy.integrate import solve_ivp

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


@nb.njit(error_model="numpy")
def _number_flux_rhs(distance, excess, thickness, above, below, flux_above,
                     particle_density, seed_mass, charging_radius, maximum_mass):
    """Derivative of F-C/m_max; C is analytic and the cap is a fixed zero bound."""
    fraction = distance/thickness
    state = (1.0 - fraction)*above + fraction*below
    mass_flux = flux_above + distance*(above[4] + .5*fraction*(below[4] - above[4]))
    floor = mass_flux/maximum_mass
    number_flux = max(excess[0], 0.0) + floor
    loss = _number_flux_loss(number_flux, mass_flux, state, particle_density, charging_radius)
    derivative = state[4]*(1.0/seed_mass - 1.0/maximum_mass) - loss
    return np.array([derivative])


def _integrate_number_flux(number_flux, thickness, state_above, state_below,
                           flux_above, particle_density, seed_mass, charging_radius,
                           rtol, atol, maximum_mass=np.inf, method='LSODA'):
    """Advance one interval with SciPy, preserving analytic mass production.

    Number flux and atol have units particles/(cm^2 s); distance is in cm.
    The integrated excess F-C/m_max is nonnegative. At a radius cap its
    entry and release are located with events. While capped, F=C/m_max
    is followed analytically; fresh production can release the constraint.
    """
    args = (thickness, state_above, state_below, flux_above,
            particle_density, seed_mass, charging_radius, maximum_mass)
    initial = max(number_flux - flux_above/maximum_mass, 0.0)
    position = 0.0

    def hits_cap(s, y, *args):
        return y[0]

    hits_cap.terminal = True
    hits_cap.direction = -1

    def release_cap(s, y, *args):
        drive = _number_flux_rhs(s, np.zeros(1), *args)[0]
        # A zero drive keeps a particle population pinned; only a transition
        # to positive growth in number flux releases the radius constraint.
        return drive if drive != 0.0 else -np.finfo(float).tiny

    release_cap.terminal = True
    release_cap.direction = 1

    def pinned(s, y, *args):
        return np.zeros(1)

    for transition in range(1000):
        if np.isfinite(maximum_mass) and initial == 0.0:
            drive = _number_flux_rhs(position, np.zeros(1), *args)[0]
            if drive <= 0.0:
                # While capped, F=C/m_max is analytic. The auxiliary solve
                # locates a release event; limit its steps within this single
                # linear-profile interval to resolve changes in the drive.
                solution = solve_ivp(
                    pinned, (position, thickness), [0.0], method=method, args=args,
                    events=release_cap, max_step=thickness/16.0,
                    rtol=rtol, atol=atol,
                )
                if not solution.success:
                    raise RuntimeError(f'{method} cap-release search failed: {solution.message}')
                position = solution.t[-1]
                if solution.status == 0:
                    break
        solution = solve_ivp(
            _number_flux_rhs, (position, thickness), [initial], method=method,
            args=args, events=hits_cap if np.isfinite(maximum_mass) else None,
            rtol=rtol, atol=atol,
        )
        if not solution.success:
            raise RuntimeError(f'{method} number-flux integration failed: {solution.message}')
        position = solution.t[-1]
        initial = max(solution.y[0, -1], 0.0)
        if solution.status == 0:
            break
        initial = 0.0
    else:
        raise RuntimeError('Too many radius-cap transitions in one interval.')
    mass_flux = flux_above + .5*(state_above[4] + state_below[4])*thickness
    flux = initial + mass_flux/maximum_mass
    return min(flux, mass_flux/seed_mass)


def _validate_haze_profile(profile, size):
    """Validate the dimensions, length, and finite values of one point profile."""
    if profile.ndim != 1 or profile.size != size:
        raise ValueError('All input profiles must be one-dimensional and match altitude.')
    if not np.all(np.isfinite(profile)):
        raise ValueError('All input profiles must contain finite values.')


def microphysics(
    altitude,
    temperature,
    air_density,
    air_number_density,
    gravity,
    production_rate,
    particle_density=1.0,
    monomer_radius=1.0e-7,
    charging_radius=9.0e-6,
    column_production=None,
    max_radius=None,
    sweep_below=None,
    rtol=1.0e-6,
    atol=1.0e-10,
    method='LSODA',
):
    """Solve steady monodisperse haze microphysics on a supplied altitude
    grid, with production of new primary particles. 
    
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
    monomer_radius : float, optional
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
        Maximum particle radius in cm. Must be finite and at least monomer_radius.
        Coagulation is limited at the imposed cap, conserving mass. Fresh
        primary-particle production can reduce radius below the cap again.
    sweep_below : float, optional
        Finite altitude in cm below which haze number density and mass flux
        are zero. Integration stops at the lowest retained sample, and a point
        exactly at the threshold is retained.
        Radius below it is held at the lowest retained sample's radius, or at
        monomer_radius if no samples remain. Represents prescribed haze removal.
    rtol : float, optional
        Relative number-flux tolerance for SciPy integration; default 1e-6.
    atol : float, optional
        Absolute number-flux tolerance in particles/(cm^2 s); default 1e-10.
    method : str, optional
        SciPy integration method; default LSODA. RK45, Radau, and BDF
        are also supported.

    Returns
    -------
    particle_radius : ndarray, shape (N,)
        Particle radius in cm, ordered from bottom to top. An empty atmosphere
        retains monomer_radius as a placeholder, except where sweep_below holds
        radius constant below the lowest retained sample.
    particle_number_density : ndarray, shape (N,)
        Particle number density in particles/cm^3, ordered from bottom to top.
    haze_mass_flux : ndarray, shape (N,)
        Positive downward mass flux in g/(cm^2 s), ordered from bottom to top.
        Zero at the upper boundary and below sweep_below. Without sweeping,
        the bottom value is total internal column production.
    """

    if method not in ('RK45', 'RK23', 'DOP853', 'Radau', 'BDF', 'LSODA'):
        raise ValueError('Unknown SciPy integration method.')

    # Check input floats
    particle_density = float(particle_density)
    monomer_radius = float(monomer_radius)
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
    if not np.isfinite(monomer_radius) or monomer_radius <= 0.0:
        raise ValueError('monomer_radius must be finite and positive.')
    if np.isnan(charging_radius) or charging_radius <= 0.0:
        raise ValueError('charging_radius must be positive (or np.inf).')
    radius_cap = np.inf
    if max_radius is not None:
        radius_cap = float(max_radius)
        if not np.isfinite(radius_cap) or radius_cap < monomer_radius:
            raise ValueError('max_radius must be finite and at least monomer_radius.')
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

    seed_mass = (4.0/3.0)*np.pi*particle_density*monomer_radius**3
    if not np.isfinite(seed_mass) or seed_mass <= 0.0:
        raise ValueError('Primary-particle mass must be finite and positive.')
    if not np.isfinite(haze_mass_flux[0]/seed_mass):
        raise ValueError('Integrated primary-particle production must be finite.')
    maximum_mass = np.inf
    if np.isfinite(radius_cap):
        maximum_mass = (4.0/3.0)*np.pi*particle_density*radius_cap**3
        if not np.isfinite(maximum_mass):
            raise ValueError('Maximum particle mass must be finite.')

    particle_radius = np.full(altitude.size, monomer_radius, dtype=np.float64)
    particle_number_density = np.zeros(altitude.size, dtype=np.float64)
    first_retained = np.searchsorted(altitude, sweep_altitude, side='left')
    if haze_mass_flux[0] == 0.0 or first_retained >= altitude.size:
        haze_mass_flux.fill(0.0)
        return particle_radius, particle_number_density, haze_mass_flux

    number_flux = 0.0
    for i in range(altitude.size - 2, first_retained - 1, -1):
        if haze_mass_flux[i] == 0.0:
            continue
        if radius_cap == monomer_radius:
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
            number_flux = _integrate_number_flux(
                number_flux, dz, state_above, state_below, haze_mass_flux[i + 1],
                particle_density, seed_mass, charging_radius, rtol, atol, maximum_mass, method
            )
        if not np.isfinite(number_flux) or number_flux <= 0.0:
            raise RuntimeError(f'Invalid particle number flux at altitude index {i}.')
        radius = np.cbrt(3.0*(haze_mass_flux[i]/number_flux)/(4.0*np.pi*particle_density))
        particle_radius[i] = max(monomer_radius, min(radius, radius_cap))
        speed = fall_velocity(gravity[i], particle_radius[i], particle_density,
                              air_density[i], air_number_density[i], temperature[i])
        particle_number_density[i] = number_flux/speed
        if not np.isfinite(particle_number_density[i]) or particle_number_density[i] <= 0.0:
            raise RuntimeError(f'Invalid particle number density at altitude index {i}.')

    if first_retained > 0:
        particle_radius[:first_retained] = particle_radius[first_retained]
        haze_mass_flux[:first_retained] = 0.0
    return particle_radius, particle_number_density, haze_mass_flux
