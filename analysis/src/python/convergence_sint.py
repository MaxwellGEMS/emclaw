import os
from convergence import Errors1D

matpath = '/simdesk/sandbox/emclaw/analysis/src/matlab/results'
matsrc  = 'analytic_centers_simpson_sint_all_hd_exact_65536'

testdir = '/simdesk/sandbox/emclaw/results/1D/convergence_sint'
basedir = '_output_'
basemin = 7
basemax = 15
frame   = 5

savedir = '/simdesk/sandbox/emclaw/results/1D/convergence_sint_test/summary_simp'

error = Errors1D(testdir,basedir,savedir,frame)

error.matsrc  = os.path.join(matpath,matsrc)
error.finesrc = os.path.join(testdir,basedir+'16')
error.basemin = basemin
error.basemax = basemax
error.debug   = True
error.p_line_range = [1,7]

error.convergence()
