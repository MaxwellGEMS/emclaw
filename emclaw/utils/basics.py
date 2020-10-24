import numpy as np
import pickle
import os

def grid_basic(dim_lims,
               cfl, co = 1.0, v = 1.0):
    """grid_basic(x_lower, x_upper, mx, cfl, co = 1, v = 1)

    Args:
        x_lower (dbl): value of x_lower (grid)
        x_upper (dbl): value of x_upper (grid)
        mx (int): number of cells
        cfl (dbl): desired CFL number
        co (dbl, optional): speed of light, from material.co. Defaults to 1.
        v (dbl, optional): speed of wave, from source.v. Defaults to 1.

    Returns:
        dx (dbl): grid space in x
        dt (dbl): time step
        tf (dbl): final time
    """

    dj = [(ds[1] - ds[0])/ds[2] for ds in dim_lims]
    dt_div = np.sum([1.0/x**2 for x in dj])
    dt = 0.9*cfl/(co*np.sqrt(dt_div))
    tf = 2.0*(dim_lims[0][1] - dim_lims[0][0])/np.min(v)

    return *dj, dt, tf

def set_outdirs(material, source, outdir, debug):
    """set_outdirs(material, source, outdir, debug)

    Args:
        material (instance): instance of material class
        source (instance): instance of source class
        outdir (string): directory where to save the data
        debug (bool): whether to do _dump_to_latex for material and source
    """
    material._outdir = outdir
    source._outdir   = outdir
    if debug: 
        material._dump_to_latex()
        source._dump_to_latex()
    return