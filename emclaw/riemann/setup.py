#!/usr/bin/env python

# How to use this file
# python setup.py build_ext -i

import os
from os.path import join as pjoin

solvers = {}
solvers['1d'] = [
    'maxwell_1d_fwave.f90',
    'maxwell_1d_nc_fwave.f90',
    'maxwell_1d_nc_tfluct.f90',
    'maxwell_1d_tfluct.f90'
]

solvers['2d'] = [
    'maxwell_2d_fwave.f90',
    'maxwell_2d_nc_fwave.f90',
    'maxwell_2d_nc_tfluct.f90',
    'maxwell_2d_tfluct.f90'
]

solvers['3d'] = [
    'maxwell_3d_nc_fwave.f90',
    'maxwell_3d_nc_tfluct.f90'
]


def configuration(package_name = 'riemann', parent_name = 'emclaw', top_path = None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration(package_name = package_name,
                           parent_name = parent_name,
                           top_path = top_path)
    this_dir = os.path.dirname(os.path.realpath(__file__))

    for dimension in range(3):
        d = str(dimension + 1) + 'd'
        for solver in solvers[d]:
            rp_subdir = 'maxwell_' + d
            ext_name = solver.split('.')[0].replace('fwave', 'rp')
            config.add_extension(name = ext_name, sources = pjoin(this_dir, rp_subdir, solver))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
