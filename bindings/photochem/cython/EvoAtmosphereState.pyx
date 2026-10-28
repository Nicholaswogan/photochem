cimport EvoAtmosphereState_pxd as state_pxd

# These private layouts match the scalar buffers in EvoAtmosphereState_wrapper.f90.
_STATE_INT_FIELDS = {
    'nq': 1, 'ng': 2, 'nsl': 3,
    'profile_enabled': 6, 'profile_mode': 7,
    'profile_hydro_pressure': 8, 'toa_enabled': 9,
}
_STATE_REAL_FIELDS = (
    'planet_mass', 'bottom_atmos', 'top_atmos', 'trop_alt', 'planet_radius',
    'profile_trop_p', 'profile_temperature_tol', 'profile_edd_tol', 'profile_extreme_factor',
    'toa_target_pressure', 'toa_pressure_factor', 'toa_extreme_pressure_factor',
)
_STATE_ARRAY_FIELDS = (
    'z', 'dz', 'temperature', 'edd', 'particle_radius', 'usol',
    'profile_pressure', 'profile_temperature', 'profile_edd',
)
_STATE_PROFILE_FIELDS = _STATE_ARRAY_FIELDS[6:]
_STATE_FLAG_FIELDS = ('profile_enabled', 'profile_hydro_pressure', 'toa_enabled')
_STATE_KEYS = frozenset(
    ('schema_version', 'species_names') + tuple(_STATE_INT_FIELDS) +
    _STATE_REAL_FIELDS + _STATE_ARRAY_FIELDS +
    tuple(key + '_allocated' for key in _STATE_PROFILE_FIELDS)
)


def _state_scalar(data, key, integer=False):
    value = data[key]
    if not isinstance(value, np.ndarray) or value.ndim != 0:
        raise PhotoException('State field "' + key + '" must be a zero-dimensional NumPy array.')
    if integer:
        if value.dtype.kind not in 'iub':
            raise PhotoException('State field "' + key + '" must contain an integer.')
        result = int(value)
        limits = np.iinfo(np.intc)
        if result < limits.min or result > limits.max:
            raise PhotoException('State field "' + key + '" exceeds the supported integer range.')
        return result
    if value.dtype.kind not in 'fiu':
        raise PhotoException('State field "' + key + '" must contain a real number.')
    return float(value)


def _state_real_array(data, key, ndim):
    value = data[key]
    if not isinstance(value, np.ndarray) or value.ndim != ndim or value.dtype.kind not in 'fiu':
        raise PhotoException('State field "' + key + '" must be a real NumPy array with ' + str(ndim) + ' dimensions.')
    return np.array(value, dtype=np.float64, order='F', copy=True)


def _state_filename(filename):
    try:
        value = os.fsdecode(os.fspath(filename))
    except TypeError as exc:
        raise PhotoException('State filename must be a string or path-like object.') from exc
    if not value or '\x00' in value or len(value.encode('utf-8')) >= 10000:
        raise PhotoException('State filename must be nonempty, contain no NUL, and be shorter than 10000 UTF-8 bytes.')
    return value


cdef class EvoAtmosphereState:
    """Independent atmospheric snapshot for restoration and storage.

    Construct from an HDF5 file or use from_dict to import a flat dictionary
    of NumPy arrays. States own their data and remain valid after the source
    model is changed or destroyed. They exclude solver history, integration
    time, boundary conditions, callbacks, and general solver settings.

    Parameters
    ----------
    filename : str or path-like
        HDF5 atmospheric snapshot to load.
    """

    cdef state_pxd.EvoAtmosphereState *_ptr
    cdef bool _ready

    def __cinit__(self, *args, **kwargs):
        self._ready = False
        self._ptr = state_pxd.allocate_evoatmospherestate()

    def __dealloc__(self):
        if self._ptr != NULL:
            state_pxd.deallocate_evoatmospherestate(self._ptr)
            self._ptr = NULL

    def __init__(self, filename):
        self.initialize_from_h5(filename)

    cdef void _require_ready(self) except *:
        if not self._ready:
            raise PhotoException('EvoAtmosphereState has not been initialized.')

    def initialize_from_h5(self, filename):
        """Load a snapshot, preserving this state on reported load errors.

        Unexpected HDF5 I/O failures can still terminate through h5fortran.

        Parameters
        ----------
        filename : str or path-like
            HDF5 atmospheric snapshot.

        Returns
        -------
        None
        """
        cdef bytes path = pystring2cstring(_state_filename(filename))
        cdef char err[ERR_LEN+1]
        state_pxd.evoatmospherestate_initialize_from_h5_wrapper(self._ptr, path, err)
        if err[0] != 0:
            raise PhotoException(err.decode('utf-8'))
        self._ready = True

    def save_to_h5(self, filename):
        """Save a snapshot through a temporary sibling file.

        Unexpected HDF5/filesystem failures can still terminate the process
        through h5fortran; predictable input errors raise PhotoException.

        Parameters
        ----------
        filename : str or path-like
            Destination file, replaced after the new file is written and closed.

        Returns
        -------
        None
        """
        self._require_ready()
        cdef bytes path = pystring2cstring(_state_filename(filename))
        cdef char err[ERR_LEN+1]
        state_pxd.evoatmospherestate_save_to_h5_wrapper(self._ptr, path, err)
        if err[0] != 0:
            raise PhotoException(err.decode('utf-8'))

    def to_dict(self):
        """Export independent copies as a flat dictionary of NumPy arrays.

        Scalars are zero-dimensional arrays; species names are fixed-width
        UTF-8 byte strings (S20), suitable for NumPy and HDF5 storage.
        Layer arrays run bottom to top. usol has shape (nq,nz) and
        particle_radius has shape (nq-ng+nsl,nz). Units follow the
        Fortran state: cm, g, K, dyn/cm², cm²/s, and molecules/cm³.
        Profile array allocation is represented by *_allocated integer
        flags and empty arrays when unallocated. Schema version is 1.

        Returns
        -------
        dict of str to numpy.ndarray
            Primary fields, compatibility metadata, and profile configuration.
        """
        self._require_ready()
        cdef ndarray integers = np.empty(10, dtype=np.intc)
        cdef ndarray reals = np.empty(12, dtype=np.float64)
        cdef ndarray lengths = np.empty(3, dtype=np.intc)
        state_pxd.evoatmospherestate_metadata_get(
            self._ptr, <int *>integers.data, <double *>reals.data, <int *>lengths.data)
        cdef int nspecies = integers[0]
        cdef ndarray names = np.empty(S_STR_LEN * nspecies, dtype=np.uint8)
        state_pxd.evoatmospherestate_species_names_get(self._ptr, &nspecies, <char *>names.data)
        raw_names = names.tobytes()
        result = {
            'schema_version': np.array(1, dtype=np.int32),
            'species_names': np.array([
                raw_names[i*S_STR_LEN:(i+1)*S_STR_LEN].rstrip(b' ')
                for i in range(nspecies)
            ], dtype='S20'),
        }
        for key, index in _STATE_INT_FIELDS.items():
            result[key] = np.array(integers[index], dtype=np.int32)
        for index, key in enumerate(_STATE_REAL_FIELDS):
            result[key] = np.array(reals[index], dtype=np.float64)
        shapes = [(int(integers[4]),)] * 4
        shapes += [(int(integers[5]), int(integers[4])), (int(integers[1]), int(integers[4]))]
        shapes += [(max(0, int(length)),) for length in lengths]
        cdef ndarray values
        cdef int field, n
        for index, (key, shape) in enumerate(zip(_STATE_ARRAY_FIELDS, shapes)):
            values = np.empty(shape, dtype=np.float64, order='F')
            field = index + 1
            n = values.size
            if n:
                state_pxd.evoatmospherestate_array_get(self._ptr, &field, &n, <double *>values.data)
            result[key] = values
        for index, key in enumerate(_STATE_PROFILE_FIELDS):
            result[key + '_allocated'] = np.array(int(lengths[index] >= 0), dtype=np.int32)
        return result

    @classmethod
    def from_dict(cls, data):
        """Construct an independent state from a schema-1 array dictionary.

        Required keys are those returned by to_dict. Missing or unknown
        fields, invalid scalar/array shapes, and unsupported versions fail
        explicitly. No defaults are supplied for incomplete saved states.

        Parameters
        ----------
        data : dict of str to numpy.ndarray
            Flat state representation, with zero-dimensional scalar arrays.

        Returns
        -------
        EvoAtmosphereState
            Validated state with independent copies of all input arrays.
        """
        cdef EvoAtmosphereState result = cls.__new__(cls)
        result._initialize_from_dict(data)
        return result

    cdef void _initialize_from_dict(self, object data) except *:
        cdef ndarray integers = np.zeros(10, dtype=np.intc)
        cdef ndarray reals = np.empty(12, dtype=np.float64)
        cdef ndarray lengths = np.empty(3, dtype=np.intc)
        cdef ndarray names, z, dz, temperature, edd, particle_radius, usol, pp, pt, pe
        cdef char err[ERR_LEN+1]

        if not isinstance(data, dict):
            raise PhotoException('Atmospheric state must be a dictionary of NumPy arrays.')
        missing = _STATE_KEYS - set(data)
        unknown = set(data) - _STATE_KEYS
        if missing or unknown:
            raise PhotoException('Invalid state dictionary keys; missing: ' +
                                 ', '.join(sorted(missing)) + '; unknown: ' +
                                 ', '.join(sorted(map(str, unknown))))
        if _state_scalar(data, 'schema_version', integer=True) != 1:
            raise PhotoException('Unsupported atmospheric state dictionary schema version.')
        for key, index in _STATE_INT_FIELDS.items():
            integers[index] = _state_scalar(data, key, integer=True)
        for index, key in enumerate(_STATE_REAL_FIELDS):
            reals[index] = _state_scalar(data, key)
        for key in _STATE_FLAG_FIELDS:
            if _state_scalar(data, key, integer=True) not in (0, 1):
                raise PhotoException('State flag "' + key + '" must be 0 or 1.')

        species = data['species_names']
        if not isinstance(species, np.ndarray) or species.ndim != 1 or species.dtype.kind not in 'SU':
            raise PhotoException('State species_names must be a one-dimensional NumPy string array.')
        encoded = []
        for value in species:
            try:
                name = bytes(value).decode('utf-8') if species.dtype.kind == 'S' else str(value)
                name_bytes = name.encode('utf-8')
            except UnicodeError as exc:
                raise PhotoException('State species_names must contain valid UTF-8 strings.') from exc
            if not name_bytes or len(name_bytes) > S_STR_LEN or b'\x00' in name_bytes:
                raise PhotoException('State species names must contain 1 to 20 UTF-8 bytes and no NUL.')
            encoded.append(name_bytes.ljust(S_STR_LEN, b' '))
        integers[0] = len(encoded)
        names = np.frombuffer(b''.join(encoded), dtype=np.uint8).copy()
        z = _state_real_array(data, 'z', 1)
        dz = _state_real_array(data, 'dz', 1)
        temperature = _state_real_array(data, 'temperature', 1)
        edd = _state_real_array(data, 'edd', 1)
        particle_radius = _state_real_array(data, 'particle_radius', 2)
        usol = _state_real_array(data, 'usol', 2)
        integers[4] = z.size
        integers[5] = particle_radius.shape[0]
        if integers[0] < 3 or integers[1] < 1 or integers[4] < 2:
            raise PhotoException('State requires species metadata, evolved gases, and at least two layers.')
        if (dz.size != z.size or temperature.size != z.size or edd.size != z.size or
                particle_radius.shape[1] != z.size or usol.shape[0] != integers[1] or usol.shape[1] != z.size):
            raise PhotoException('State atmospheric array shapes are inconsistent.')
        pp = _state_real_array(data, 'profile_pressure', 1)
        pt = _state_real_array(data, 'profile_temperature', 1)
        pe = _state_real_array(data, 'profile_edd', 1)
        for index, (key, values) in enumerate(zip(_STATE_PROFILE_FIELDS, (pp, pt, pe))):
            allocated = _state_scalar(data, key + '_allocated', integer=True)
            if allocated not in (0, 1):
                raise PhotoException('State allocation flags must be 0 or 1.')
            if not allocated and values.size:
                raise PhotoException('Unallocated state profile arrays must be empty: ' + key)
            lengths[index] = values.size if allocated else -1

        state_pxd.evoatmospherestate_create_wrapper(
            self._ptr, <int *>integers.data, <double *>reals.data, <int *>lengths.data,
            <char *>names.data, <double *>z.data, <double *>dz.data,
            <double *>temperature.data, <double *>edd.data, <double *>particle_radius.data,
            <double *>usol.data, <double *>pp.data, <double *>pt.data, <double *>pe.data, err)
        if err[0] != 0:
            raise PhotoException(err.decode('utf-8'))
        self._ready = True

