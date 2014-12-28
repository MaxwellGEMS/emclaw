import os
from convergence import Errors1D

matpath = '../matlab/results'
matsrc  = '_sin_nc_65536.mat'

testdir = '/simdesk/sandbox/emclaw/results/1D/_convergence_sin_homogeneous_tfluct'
compare_dir = '/simdesk/sandbox/emclaw/results/1D/_convergence_sin'
basedir = '_output_'
basemin = 7
basemax = 15
frame   = 5

savedir = os.path.join(testdir,'_summaryb')

error = Errors1D(testdir,basedir,savedir,frame)

error.matsrc  = os.path.join(matpath,matsrc)
error.finesrc = os.path.join(testdir,basedir+'16')
error.basemin = basemin
error.basemax = basemax
error.debug   = True
error.p_line_range = [1,7]
error.homogeneous = True
error.compare_dir = compare_dir
error.compare = True

error.convergence()
