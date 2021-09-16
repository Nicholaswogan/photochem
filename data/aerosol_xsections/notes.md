
# Aerosol optical properties


**Folder and File formats**

Each folder corresponds to optical properties of different types of aerosols. For example `khare1984` contains the measured index of refraction of laboratory synthetic hydrocarbon aerosols analogous to the ones that are made in Titan's atmosphere.

Each folder should contain a few files: `mie.dat` (required), `frac.dat` (optional), and an ascii file containing measured index of refractions.

`mie.dat` - Fortran binary file containing calculated optical properties using mie theory.

`frac.dat` (optional) - Fortran binary file containing calculated optical properties using fractal theory following [Wolf and Toon (2010)](https://science.sciencemag.org/content/328/5983/1266.abstract).

## Folders

`khare1984` is data for optical properties of fractal hydrocarbon aerosols from [Khare et al. (1984)](https://www.sciencedirect.com/science/article/pii/0019103584901428) downloaded from hitran.org at this link: https://hitran.org/data/Aerosols-2016/ascii/exoplanets/khare_tholins.dat.

`palmer1975` is data for optical properties of H2SO4 from [Palmer et al. (1975)](https://www.osapublishing.org/ao/abstract.cfm?uri=ao-14-1-208) downloaded from hitran.org at this link: https://hitran.org/data/Aerosols-2016/ascii/single_files/palmer_williams_h2so4.dat. I used the column with 95.6% H2SO4. 






