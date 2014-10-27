import os
from convergence import Errors1D

matpath = '../matlab/results'
matsrc  = '_rip_nc_65536.mat'

testdir = '/simdesk/sandbox/emclaw/results/1D/_convergence_rip'
basedir = '_output_'
basemin = 7
basemax = 15
frame   = 61

savedir = '/simdesk/sandbox/emclaw/results/1D/_convergence_rip/_summary'

error = Errors1D(testdir,basedir,savedir,frame)

error.matsrc  = os.path.join(matpath,matsrc)
error.finesrc = os.path.join(testdir,basedir+'16')
error.basemin = basemin
error.basemax = basemax
error.debug   = True
error.p_line_range = [2,7]

error.convergence()
