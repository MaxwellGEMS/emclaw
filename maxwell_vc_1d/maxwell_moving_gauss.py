import sys
import os
import numpy as np

sys.path.append(os.path.realpath('../utils'))
sys.path.append(os.path.realpath('../'))


from utils.materials import Material1D
from utils.sources import Source1D

material = Material1D()
material.bkg_er = 1.5
material.bkg_mr = 1.5
material.shape = 'moving_gauss'
material.setup()

source = Source1D(material,shape='pulse',wavelength=2.0)
source.offset = -5.0
source.setup()

x_lower = 0.
x_upper = 400.

def grid_basic(x_lower,x_upper,mx,cfl):
    dx = (x_upper-x_lower)/mx
    dt = 0.9*cfl/(material.co*np.sqrt(1.0/(dx**2)))
    tf = (x_upper-x_lower)/source.v

    return dx,dt,tf

def em1D(mx=2**11,num_frames=100,cfl=1.0,outdir='./_output',before_step=True,debug=False,velocity=0.59):
    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    if MPI.COMM_WORLD.rank==0:
        material.dump()
        source.dump()

    material.velocity_e = velocity
    material.velocity_m = velocity

    num_eqn   = 2
    num_waves = 2
    num_aux   = 4

#   grid pre calculations and domain setup
    dx,dt,tf = grid_basic(x_lower,x_upper,mx,cfl)
    x = pyclaw.Dimension('x',x_lower,x_upper,mx)
    domain = pyclaw.Domain([x])

#   Solver settings
    solver=pyclaw.SharpClawSolver1D()
    solver.num_waves  = num_waves
    solver.num_eqn    = num_eqn
    solver.weno_order = 5

    solver.dt_variable = True
    solver.dt_initial  = dt/2.0
    solver.dt_max      = dt
    solver.max_steps   = int(2*tf/dt)
    
#   Import Riemann and Tfluct solvers
    import maxwell_1d_rp
    import maxwell_1d_tfluct

    solver.tfluct_solver = True
    solver.fwave = True
    
    solver.rp = maxwell_1d_rp
    solver.tfluct = maxwell_1d_tfluct
    
    solver.cfl_max = cfl+0.5
    solver.cfl_desired = cfl

#   boundary conditions
    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.user_bc_lower = source.scattering_bc

    solver.aux_bc_lower[0]= pyclaw.BC.custom
    solver.aux_bc_upper[0]= pyclaw.BC.custom
    solver.user_aux_bc_lower = material.setaux_lower
    solver.user_aux_bc_upper = material.setaux_upper

#   before step configure
    if before_step:
            solver.call_before_step_each_stage = True
            solver.before_step = material.update_aux

#   state setup
    state = pyclaw.State(domain,num_eqn,num_aux)
    
    state.problem_data['chi2_e'] = material.chi2_e
    state.problem_data['chi3_e'] = material.chi3_e
    state.problem_data['chi2_m'] = material.chi2_m
    state.problem_data['chi3_m'] = material.chi3_m
    state.problem_data['eo'] = material.eo
    state.problem_data['mo'] = material.mo
    state.problem_data['co'] = material.co
    state.problem_data['zo'] = material.zo
    state.problem_data['dx'] = state.grid.x.delta

#   array initialization
    source.init(state)
    material.init(state)

#   controller
    claw = pyclaw.Controller()
    claw.tfinal = tf
    claw.num_output_times = num_frames
    claw.solver = solver
    claw.solution = pyclaw.Solution(state,domain)
    claw.outdir = outdir
    claw.write_aux_always = True

    return claw

if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em1D)
