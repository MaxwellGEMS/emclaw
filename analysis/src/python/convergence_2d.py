import os
from convergence import Errors2D

matpath = '/media/noor/simdesk/results/analysis/results'
matsrc  = 'analytic_centers_simpson__all_hd_exact_65536.mat'

testdir = '/media/shaheen/dev/results/maxwell_2d/convergence_averaged'
basedir = '_output_plane_'
basemin = 7
basemax = 13
frame   = 45

savedir = '/simdesk/sandbox/emclaw/results/2D/convergence_averaged_test/summary'

error = Errors2D(testdir,basedir,savedir,frame)

error.finesrc = os.path.join(testdir,basedir+'13')
error.basemin = basemin
error.basemax = basemax
error.debug   = True

error.convergence()
