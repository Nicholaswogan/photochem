from ._clima import AdiabatClimate, ClimaException
from ._clima import (
    RCE_SOLVE_HYBRJ_ONLY,
    RCE_SOLVE_PTC_THEN_HYBRJ,
    RCE_SOLVE_HYBRJ_THEN_PTC_THEN_HYBRJ,
)
from ._clima import rebin, rebin_with_errors

# Limits OpenMP threads to 1
from threadpoolctl import threadpool_limits
_ = threadpool_limits(limits=1, user_api='openmp')
