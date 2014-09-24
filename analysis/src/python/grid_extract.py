from scipy.io import loadmat,savemat
from glob import glob
import os
from clawpack.pyclaw import Solution


testpath = '/simdesk/sandbox/emclaw/maxwell_vc_1d'
testbase = '_gauss_v059_pulse_'
testdirs = sorted(glob(os.path.join(testpath,testbase+'*')))


claw_grid = {}


for m, dirs in enumerate(testdirs):
    sol1 = Solution()

    # load the pyclaw solution, q
    sol1.read(75,path=dirs,file_format='petsc',read_aux=False)

    claw_grid['grid'+str(m)] = {}
    claw_grid['grid'+str(m)]['centers'] = sol1.state.grid.x.centers
    claw_grid['grid'+str(m)]['edges'] = sol1.state.grid.x.edges
    claw_grid['delta'+str(m)] = sol1.state.grid.x.delta


claw_grid['t'] = sol1.state.t
savemat('claw_grids75',claw_grid)

