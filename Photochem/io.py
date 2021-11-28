from scipy.io import FortranFile, FortranEOFError
import numpy as np

def read_evolve_output(filename):
    f = FortranFile(filename, 'r')
    # metadata
    nq = f.read_record(np.int32)[0]
    nz = f.read_record(np.int32)[0]
    z = f.read_record(np.double)
    tmp = f.read_record(np.dtype("S20"))
    species_names = np.char.strip(tmp.astype('U')).tolist()
    nt = f.read_record(np.int32)[0]

    # results
    time = np.empty((nt),order='F')
    usol = np.empty((nq,nz,nt),order='F')
    err = False
    for i in range(nt):
        try:
            time[i] = f.read_record(np.double)[0]
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
        usol = np.delete(usol,slice(i,nt),axis=2)

    # return dictionary
    sol = {}
    sol['species_names'] = species_names
    sol['alt'] = z/1e5
    sol['time'] = time
    sol['usol'] = usol
    return sol
    
def reformat_output_dict(sol):
    out = {}
    out['time'] = sol['time']
    out['alt'] = sol['alt']
    for i,sp in enumerate(sol['species_names']):
        out[sp] = sol['usol'][i,:,:]
    return out
    
