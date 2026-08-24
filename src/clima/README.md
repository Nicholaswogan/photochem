# Imported Clima source

The sources in this directory were imported from Clima v0.7.5 at commit
`e37e32c0904c4651faac440a63652fca9b17622f`. Clima is distributed under the
GNU General Public License v3.0, the same license used by this repository.

The files formerly stored under Clima's `src/adiabat/` and `src/radtran/`
directories are flattened here. The experimental `src/climate/` implementation
was deliberately excluded and the public Fortran facade therefore does not
export its `Climate` type.

The vendored `linear_interpolation_module.F90` dependency is kept under
`src/dependencies/` with its original BSD license notice and attribution.
