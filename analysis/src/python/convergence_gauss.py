import os
from convergence import Errors1D

matpath = '/media/noor/simdesk/results/analysis/results'
matsrc  = 'analytic_centers_simpson__all_hd_exact_65536.mat'

testdir = '/simdesk/sandbox/emclaw/results/1D/_convergence_alt_src_averaged'
basedir = '_output_'
basemin = 7
basemax = 15
frame   = 61

savedir = '/simdesk/sandbox/emclaw/results/1D/convergence_gauss_test/summary_simp'

error = Errors1D(testdir,basedir,savedir,frame)

error.matsrc  = os.path.join(matpath,matsrc)
error.finesrc = os.path.join(testdir,basedir+'16')
error.basemin = basemin
error.basemax = basemax
error.debug   = True
error.p_line_range = [2,7]

error.convergence()
