# Utilities API

These helpers prepare common input files and stellar spectra or read model output. Functions that download external data require network access and are identified in their docstrings.

## Climate input files

::: photochem.utils.climate
    options:
      members:
        - species_dict_for_climate
        - species_file_for_climate
        - settings_dict_for_climate
        - settings_file_for_climate
      show_root_heading: false

## Stellar spectra

::: photochem.utils.stars
    options:
      members:
        - stefan_boltzmann
        - stefan_boltzmann_inverse
        - equilibrium_temperature
        - equilibrium_temperature_inverse
        - blackbody_cgs
        - blackbody
        - grid_at_resolution
        - make_bins
        - energy_in_spectrum
        - scale_spectrum_to_planet
        - append_blackbody_to_stellar_spectrum
        - photochem_spectrum_string
        - save_photochem_spectrum
        - rebin_to_needed_resolution
        - solar_spectrum
        - hazmat_spectrum
        - print_hazmat_stars
        - muscles_spectrum
        - closest_muscles_to_Teff
        - print_muscles_stars
      show_root_heading: false

## Model output

::: photochem.io
    options:
      members:
        - reformat_output_dict
        - evo_read_evolve_output
      show_root_heading: false

## Reaction mechanisms

::: photochem.utils
    options:
      members:
        - FormatReactions
        - FormatSettings
        - resave_mechanism_with_atoms
        - generate_zahnle_earth_thermo
        - zahnle_rx_and_thermo_files
        - atmos2yaml
        - atmosbc2yaml
        - vulcan2yaml
        - photochem2cantera
      show_root_heading: false
