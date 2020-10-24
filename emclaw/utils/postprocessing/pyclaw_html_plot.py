import os
import sys
import shutil
import errno
from glob import glob
from clawpack.petclaw import plot
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams.update({'font.size': 10})

def copy(src, dest):
    try:
        shutil.copytree(src, dest)
    except OSError as e:
        # If the error was caused because the source wasn't a directory
        if e.errno == errno.ENOTDIR:
            shutil.copy(src, dst)
        else:
            print('Directory not copied. Error: %s' % e)

def html_plot(src):
    plot.html_plot(outdir=src)

def main_plot(outdir='./_output',multiple=False,overwrite=False):
    if multiple:
        outdir = outdir+'*'
    
    outdirs = sorted(glob(outdir))
    for dirs in outdirs:
        if overwrite or not os.path.exists(os.path.join(outdir,'_plots')):
            plot.html_plot(outdir=dirs)
            copy('./_plots',os.path.join(dirs,'_plots'))
            shutil.rmtree('./_plots')

if __name__ == "__main__":
    from clawpack.pyclaw import util
    args,app_args = util._info_from_argv(sys.argv)
    main_plot(**app_args)

    
