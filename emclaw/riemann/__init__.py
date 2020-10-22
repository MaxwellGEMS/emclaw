#!/usr/bin/env python
# encoding: utf-8
"""
Wave propagation Riemann solvers implemented in Fortran.
"""
__pdoc__ = {
    'maxwell_1d':False,
    'maxwell_2d':False,
    'maxwell_3d':False,
    'build': False,
    '__pycache__':False
}

from . import maxwell_1d_nc_rp
from . import maxwell_1d_nc_tfluct
from . import maxwell_1d_rp
from . import maxwell_1d_tfluct

from . import maxwell_2d_nc_rp
from . import maxwell_2d_nc_tfluct
from . import maxwell_2d_rp
from . import maxwell_2d_tfluct

from . import maxwell_3d_nl_rp as maxwell_3d_nc_rp
from . import maxwell_3d_nl_tfluct as maxwell_3d_nc_tfluct
