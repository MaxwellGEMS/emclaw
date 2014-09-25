#!/usr/bin/env python

# How to use this file
# python setup.py build_ext -i

import os
from os.path import join as pjoin

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration
    config = Configuration('', parent_package, top_path)
    this_dir = os.path.dirname(os.path.realpath(__file__))

    config.add_extension('maxwell_3d_rp',pjoin(this_dir,'maxwell_3d_nl_fwave.f90'))
    config.add_extension('maxwell_3d_tfluct',pjoin(this_dir,'maxwell_3d_nl_tfluct.f90'))

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
