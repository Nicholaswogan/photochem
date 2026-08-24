from ._clima import AdiabatClimate, ClimaException
from ._clima import rebin # rebin routine from futils

# Limits OpenMP threads to 1
from threadpoolctl import threadpool_limits
_ = threadpool_limits(limits=1, user_api='openmp')
