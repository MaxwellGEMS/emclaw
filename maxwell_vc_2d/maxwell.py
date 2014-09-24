import sys
import os
import numpy as np

sys.path.append(os.path.realpath('../utils'))
sys.path.append(os.path.realpath('../'))


from utils.materials import Material2D
from utils.sources import Source2D

x_lower = 0.0
x_upper = 10.0

y_lower = 0.0
y_upper = 10.0

material = Material2D(shape='homogeneous',metal=False)
material.setup()
material._calculate_n()

source = Source2D(material,shape='pulse',wavelength=2.0)
source.offset[1] = -5.0
source.setup()

def grid_basic(x_lower,x_upper,y_lower,y_upper,mx,my,cfl):
    dx = (x_upper-x_lower)/mx
    dy = (y_upper-y_lower)/my
    dt = 0.90/(material.co*np.sqrt(1.0/(dx**2)+1.0/(dy**2)))
    tf = (x_upper-x_lower+5.0)/material.v
    
    return dx,dy,dt,tf

def em2D(mx=2**11,my=2**6,num_frames=100,cfl=1.0,outdir='./_output',before_step=True,debug=False):
    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    if MPI.COMM_WORLD.rank==0:
        material.dump()
        source.dump()

    num_eqn   = 3
    num_waves = 2
    num_aux   = 6

#   grid pre calculations and domain setup
    dx,dy,dt,tf = grid_basic(x_lower,x_upper,y_lower,y_upper,mx,my,cfl)

    x = pyclaw.Dimension('x',x_lower,x_upper,mx)
    y = pyclaw.Dimension('y',y_lower,y_upper,my)
    
    domain = pyclaw.Domain([x,y])

#   Solver settings
    solver = pyclaw.SharpClawSolver2D()
    solver.num_waves  = num_waves
    solver.num_eqn    = num_eqn
    solver.weno_order = 5

    solver.dt_variable = True
    solver.dt_initial  = dt/2.0
    solver.dt_max      = dt
    solver.max_steps   = int(2*tf/dt)
    
#   Import Riemann and Tfluct solvers
    import maxwell_2d_rp
    import maxwell_2d_tfluct

    solver.tfluct_solver = True
    solver.fwave = True
    
    solver.rp = maxwell_2d_rp
    solver.tfluct = maxwell_2d_tfluct
    
    solver.cfl_max = cfl+0.5
    solver.cfl_desired = cfl

#   boundary conditions
    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap
    solver.user_bc_lower = source.scattering_bc

    solver.aux_bc_lower[0]= pyclaw.BC.custom
    solver.aux_bc_lower[1]= pyclaw.BC.extrap
    solver.aux_bc_upper[0]= pyclaw.BC.extrap
    solver.aux_bc_upper[1]= pyclaw.BC.extrap
    
    solver.user_aux_bc_lower = material.setaux_lower
    solver.user_aux_bc_upper = material.setaux_upper

#   before step configure
    if before_step:
            solver.call_before_step_each_stage = True
            solver.before_step = material.update_aux

#   state setup
    state = pyclaw.State(domain,num_eqn,num_aux)
    
    state.problem_data['chi2']  = material.chi2
    state.problem_data['chi3']  = material.chi3
    state.problem_data['vac1']  = material.eo
    state.problem_data['vac2']  = material.eo
    state.problem_data['vac3']  = material.mo
    state.problem_data['co'] = material.co
    state.problem_data['zo'] = material.zo
    state.problem_data['dx'] = state.grid.x.delta
    state.problem_data['dy'] = state.grid.y.delta

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
    output = run_app_from_main(em2D)
