
# Aerosol optical properties

I use the `miepython` python package to compute aerosol optical properties: https://github.com/scottprahl/miepython

**Folder and File formats**

Each folder corresponds to optical properties of different types of aerosols. For example `khare1984` contains the measured index of refraction of laboratory synthetic hydrocarbon aerosols analogous to the ones that are made in Titan's atmosphere.

Each folder should contain a few files: `mie_<folder-name>.dat` (required), `frac_<folder-name>.dat` (optional), and an ascii file containing measured index of refractions.

`mie_<folder-name>.dat` - Fortran binary file containing calculated optical properties using mie theory. The records in the mie binary file have the following organization
1. `nw`, number of wavelength bins (integer(4))
2. wavelengths in nm (real(8), dimension(nw))
3. `nrad`, number of radii bins (integer(4))
4. radii in micrometer (real(8), dimension(nrad))
5. single scattering albedo (real(8), dimension(nrad, nw))
6. extinction (real(8), dimension(nrad, nw))
7. asymmetry factor (real(8), dimension(nrad, nw))

`frac_<folder-name>.dat` (optional) - Fortran binary file containing calculated optical properties using fractal theory following [Wolf and Toon (2010)](https://science.sciencemag.org/content/328/5983/1266.abstract).

## Folders

`khare1984` is data for optical properties of fractal hydrocarbon aerosols from [Khare et al. (1984)](https://www.sciencedirect.com/science/article/pii/0019103584901428) downloaded from hitran.org at this link: https://hitran.org/data/Aerosols-2016/ascii/exoplanets/khare_tholins.dat.

`palmer1975` is data for optical properties of H2SO4 from [Palmer et al. (1975)](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-14-1-208) downloaded from hitran.org at this link: https://hitran.org/data/Aerosols-2016/ascii/single_files/palmer_williams_h2so4.dat. I used the column with 95.6% H2SO4. 






