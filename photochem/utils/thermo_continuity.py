"""Inspect and repair gas thermodynamic fits in Photochem and Clima YAML data."""

from copy import deepcopy
from math import isfinite, log


RGAS = 8.31446261815324  # J / (mol K), matching photochem_const.f90


def _enthalpy_entropy(model, coefficients, temperature):
    """Match the formulas used by Photochem's Fortran Gibbs evaluators."""
    c = coefficients
    t = temperature
    if model == "Shomate":
        x = t / 1000.0
        enthalpy = 1000.0 * (
            c[0] * x + c[1] * x**2 / 2 + c[2] * x**3 / 3
            + c[3] * x**4 / 4 - c[4] / x + c[5]
        )
        entropy = (
            c[0] * log(x) + c[1] * x + c[2] * x**2 / 2
            + c[3] * x**3 / 3 - c[4] / (2 * x**2) + c[6]
        )
    elif model == "NASA7":
        enthalpy = RGAS * t * (
            c[0] + c[1] * t / 2 + c[2] * t**2 / 3
            + c[3] * t**3 / 4 + c[4] * t**4 / 5 + c[5] / t
        )
        entropy = RGAS * (
            c[0] * log(t) + c[1] * t + c[2] * t**2 / 2
            + c[3] * t**3 / 3 + c[4] * t**4 / 4 + c[6]
        )
    elif model == "NASA9":
        enthalpy = RGAS * t * (
            -c[0] / t**2 + c[1] * log(t) / t + c[2]
            + c[3] * t / 2 + c[4] * t**2 / 3
            + c[5] * t**3 / 4 + c[6] * t**4 / 5 + c[7] / t
        )
        entropy = RGAS * (
            -c[0] / (2 * t**2) - c[1] / t + c[2] * log(t)
            + c[3] * t + c[4] * t**2 / 2
            + c[5] * t**3 / 3 + c[6] * t**4 / 4 + c[8]
        )
    else:
        raise ValueError(f"Unsupported thermodynamic model: {model!r}")
    return enthalpy, entropy


def _heat_capacity(model, coefficients, temperature):
    """Evaluate Cp in J/(mol K) using the same polynomial conventions."""
    c = coefficients
    t = temperature
    if model == "Shomate":
        x = t / 1000.0
        return c[0] + c[1] * x + c[2] * x**2 + c[3] * x**3 + c[4] / x**2
    if model == "NASA7":
        return RGAS * (c[0] + c[1] * t + c[2] * t**2 + c[3] * t**3 + c[4] * t**4)
    if model == "NASA9":
        return RGAS * (c[0] / t**2 + c[1] / t + c[2] + c[3] * t
                       + c[4] * t**2 + c[5] * t**3 + c[6] * t**4)
    raise ValueError(f"Unsupported thermodynamic model: {model!r}")


def _validated_thermo(species):
    thermo = species.get("thermo")
    if thermo is None:
        return None
    model = thermo.get("model")
    expected = {"Shomate": 7, "NASA7": 7, "NASA9": 9}
    if model not in expected:
        raise ValueError(f"{species['name']}: unsupported thermodynamic model {model!r}")
    edges = thermo.get("temperature-ranges", [])
    rows = thermo.get("data", [])
    if len(edges) != len(rows) + 1 or not rows:
        raise ValueError(f"{species['name']}: temperature ranges and data do not match")
    if any(not isfinite(float(x)) for x in edges) or any(
        float(right) <= float(left) for left, right in zip(edges, edges[1:])
    ):
        raise ValueError(f"{species['name']}: temperature ranges must increase")
    if any(len(row) != expected[model] or any(not isfinite(float(x)) for x in row)
           for row in rows):
        raise ValueError(f"{species['name']}: invalid {model} coefficients")
    if any(float(x) <= 0 for x in edges[1:-1]):
        raise ValueError(f"{species['name']}: joins must have positive temperature")
    return thermo


def check_thermo_continuity(network, *, enthalpy_tolerance=1e-3,
                            entropy_tolerance=1e-6, heat_capacity_tolerance=1e-6,
                            gibbs_tolerance=1e-3):
    """Report discontinuities in gas thermodynamic polynomial joins.

    Accepts a parsed ``EvoAtmosphere`` reaction file, ``AdiabatClimate``
    species file, or ``Equilibrate`` thermodynamic input. Only entries in
    ``network["species"]`` with a ``thermo`` block are checked. Species marked
    ``condensate: true`` are skipped. Reactions, saturation data, and particles
    are ignored.

    Parameters
    ----------
    network : dict
        YAML mapping from an EvoAtmosphere reaction file, AdiabatClimate
        species file, or Equilibrate thermodynamic input.
    enthalpy_tolerance : float, optional
        Maximum allowed enthalpy jump in J/mol. Default is ``1e-3``.
    entropy_tolerance : float, optional
        Maximum allowed entropy jump in J/(mol K). Default is ``1e-6``.
    heat_capacity_tolerance : float, optional
        Maximum allowed heat-capacity jump in J/(mol K). Default is ``1e-6``.
    gibbs_tolerance : float, optional
        Maximum allowed Gibbs-energy jump in J/mol. Default is ``1e-3``.

    Returns
    -------
    list of dict
        One entry per join exceeding any tolerance. Each entry identifies the
        species, model, boundary temperature in K, zero-based segment indices,
        and right-minus-left jumps in enthalpy, entropy, heat capacity, and
        Gibbs energy. An empty list means no jumps exceed the tolerances.
    """
    discontinuities = []
    for species in network["species"]:
        if species.get("condensate") is True:
            continue
        thermo = _validated_thermo(species)
        if thermo is None:
            continue
        model = thermo["model"]
        for index, temperature in enumerate(thermo["temperature-ranges"][1:-1], 1):
            left_h, left_s = _enthalpy_entropy(model, thermo["data"][index - 1], temperature)
            right_h, right_s = _enthalpy_entropy(model, thermo["data"][index], temperature)
            delta_h = right_h - left_h
            delta_s = right_s - left_s
            delta_g = delta_h - temperature * delta_s
            delta_cp = (_heat_capacity(model, thermo["data"][index], temperature)
                        - _heat_capacity(model, thermo["data"][index - 1], temperature))
            if (abs(delta_h) > enthalpy_tolerance or
                    abs(delta_s) > entropy_tolerance or
                    abs(delta_cp) > heat_capacity_tolerance or
                    abs(delta_g) > gibbs_tolerance):
                discontinuities.append({
                    "species": species["name"],
                    "model": model,
                    "temperature": float(temperature),
                    "left_segment": index - 1,
                    "right_segment": index,
                    "delta_enthalpy": delta_h,
                    "delta_entropy": delta_s,
                    "delta_heat_capacity": delta_cp,
                    "delta_gibbs": delta_g,
                })
    return discontinuities


def make_thermo_continuous(network, *, reference_temperature=298.15):
    """Return a copy with continuous gas thermodynamic polynomial joins.

    Accepts a parsed ``EvoAtmosphere`` reaction file, ``AdiabatClimate``
    species file, or ``Equilibrate`` thermodynamic input. Species marked
    ``condensate: true`` are skipped. The segment containing
    ``reference_temperature`` is kept unchanged, and adjacent segments are
    aligned outward from it. The repair matches heat capacity, enthalpy,
    entropy, and therefore Gibbs energy at each join. It can change one
    heat-capacity coefficient and two integration constants throughout each
    adjusted segment. Reactions, saturation data, particles, and the input
    mapping are left unchanged.

    Parameters
    ----------
    network : dict
        YAML mapping from an EvoAtmosphere reaction file, AdiabatClimate
        species file, or Equilibrate thermodynamic input. Each multi-segment
        gas fit must span ``reference_temperature``.
    reference_temperature : float, optional
        Temperature in K selecting the unchanged segment for each species.
        Default is ``298.15``.

    Returns
    -------
    repaired : dict
        Deep copy of ``network`` with adjusted gas thermodynamic coefficients.
    corrections : list of dict
        One entry per adjusted segment, identifying its species, model, join
        temperature, zero-based segment index, heat-capacity shift, enthalpy
        and entropy integration-constant shifts, and the original Gibbs jump.
        Review large corrections against the source thermodynamic data.
    """
    repaired = deepcopy(network)
    corrections = []
    for species in repaired["species"]:
        if species.get("condensate") is True:
            continue
        thermo = _validated_thermo(species)
        if thermo is None or len(thermo["data"]) == 1:
            continue
        # Ensure independently editable coefficients even if input species
        # happened to share a thermodynamic dictionary object.
        thermo = deepcopy(thermo)
        species["thermo"] = thermo
        edges = thermo["temperature-ranges"]
        rows = thermo["data"]
        model = thermo["model"]
        if not edges[0] <= reference_temperature <= edges[-1]:
            raise ValueError(
                f"{species['name']}: reference temperature is outside the fit range"
            )
        anchor = next(
            (i for i in range(len(rows))
             if edges[i] <= reference_temperature < edges[i + 1]),
            len(rows) - 1,
        )
        for target in list(range(anchor + 1, len(rows))) + list(range(anchor - 1, -1, -1)):
            neighbor = target - 1 if target > anchor else target + 1
            boundary = edges[target] if target > anchor else edges[target + 1]
            neighbor_h, neighbor_s = _enthalpy_entropy(model, rows[neighbor], boundary)
            original_h, original_s = _enthalpy_entropy(model, rows[target], boundary)
            shift_cp = (_heat_capacity(model, rows[neighbor], boundary)
                        - _heat_capacity(model, rows[target], boundary))
            if (abs(shift_cp) <= 1e-12 and abs(neighbor_h - original_h) <= 1e-9
                    and abs(neighbor_s - original_s) <= 1e-12):
                continue
            if model == "Shomate":
                rows[target][0] += shift_cp
            elif model == "NASA7":
                rows[target][0] += shift_cp / RGAS
            else:  # NASA9
                rows[target][2] += shift_cp / RGAS
            target_h, target_s = _enthalpy_entropy(model, rows[target], boundary)
            shift_h = neighbor_h - target_h
            shift_s = neighbor_s - target_s
            if model == "Shomate":
                rows[target][5] += shift_h / 1000.0
                rows[target][6] += shift_s
            elif model == "NASA7":
                rows[target][5] += shift_h / RGAS
                rows[target][6] += shift_s / RGAS
            else:  # NASA9
                rows[target][7] += shift_h / RGAS
                rows[target][8] += shift_s / RGAS
            corrections.append({
                "species": species["name"],
                "model": model,
                "temperature": float(boundary),
                "adjusted_segment": target,
                "heat_capacity_shift": shift_cp,
                "enthalpy_shift": shift_h,
                "entropy_shift": shift_s,
                "gibbs_jump_before": (
                    (original_h - neighbor_h) - boundary * (original_s - neighbor_s)
                ) * (1 if target > anchor else -1),
            })
    return repaired, corrections
