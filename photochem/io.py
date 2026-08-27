from scipy.io import FortranFile, FortranEOFError
import numpy as np
   
def reformat_output_dict(sol):
    """Convert raw evolution output to a species-keyed dictionary.

    Parameters
    ----------
    sol : dict
        Output from [evo_read_evolve_output][photochem.io.evo_read_evolve_output].

    Returns
    -------
    dict[str, ndarray]
        ``time`` in seconds, ``alt`` in km, and one ``(nz, nt)`` number-density
        array in molecules/cm^3 for each evolved species.
    """
    out = {}
    out['time'] = sol['time']
    out['alt'] = sol['alt']
    for i,sp in enumerate(sol['species_names']):
        out[sp] = sol['usol'][i,:,:]
    return out
    
def evo_read_evolve_output(filename):
    """Read the unformatted binary file written by ``EvoAtmosphere.evolve``.

    A partially written final record is discarded, so output from an
    interrupted integration can still be inspected.

    Parameters
    ----------
    filename : str or path-like
        Evolution output file.

    Returns
    -------
    dict
        Species names; ``time`` in seconds; ``top_atmos`` in cm; ``alt`` in km
        with shape ``(nz, nt)``; and evolved gas and condensed-material number
        densities in molecules/cm^3 with shape ``(nq, nz, nt)``.
    """
    f = FortranFile(filename, 'r')
    # metadata
    nq = f.read_record(np.int32)[0]
    nz = f.read_record(np.int32)[0]
    tmp = f.read_record(np.dtype("S20"))
    species_names = np.char.strip(tmp.astype('U')).tolist()
    nt = f.read_record(np.int32)[0]

    # results
    time = np.empty((nt),order='F')
    top_atmos = np.empty((nt),order='F')
    z = np.empty((nz,nt),order='F')
    usol = np.empty((nq,nz,nt),order='F')
    err = False
    for i in range(nt):
        try:
            time[i] = f.read_record(np.double)[0]
        except FortranEOFError:
            err = True
            break
        try:
            top_atmos[i] = f.read_record(np.double)[0]
        except FortranEOFError:
            err = True
            break
        try:
            z[:,i] = f.read_record(np.double)
        except FortranEOFError:
            err = True
            break
        try:
            usol[:,:,i] = f.read_record(np.double).reshape((nq,nz),order='F')
        except FortranEOFError:
            err = True
            break
    f.close()

    if err:
        time = np.delete(time,slice(i,nt),axis=0)
        top_atmos = np.delete(top_atmos,slice(i,nt),axis=0)
        z = np.delete(z,slice(i,nt),axis=1)
        usol = np.delete(usol,slice(i,nt),axis=2)

    # return dictionary
    sol = {}
    sol['species_names'] = species_names
    sol['time'] = time
    sol['top_atmos'] = top_atmos
    sol['alt'] = z/1e5
    sol['usol'] = usol
    return sol
    
