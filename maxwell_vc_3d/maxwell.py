import sys
import os
import numpy as np

sys.path.append(os.path.realpath('../utils'))
sys.path.append(os.path.realpath('../'))


from utils.materials import Material3D
from utils.sources import Source3D

x_lower = 0.0
x_upper = 10.0

y_lower = 0.0
y_upper = 10.0

z_lower = 0.0
z_upper = 10.0

sy = y_upper-y_lower
sx = x_upper-x_lower
sz = z_upper-z_lower

mid_point_y = (y_upper-y_lower)/2.0
mid_point_z = (z_upper-z_lower)/2.0

material = Material3D(shape='homogeneous',metal=False)
material.setup()
material._calculate_n()

source = Source3D(material,shape='pulse',wavelength=2.0)
source.offset[1] = -5.0
source.offset[5] = -5.0
source.transversal_shape = 'bessel'
source.setup()
source.transversal_offset = [sy/2.0,sz/2.0]
source.transversal_width = [sy/2.0,sz/2.0]

def grid_basic(x_lower,x_upper,y_lower,y_upper,z_lower,z_upper,mx,my,mz,cfl):
    dx = (x_upper-x_lower)/mx
    dy = (y_upper-y_lower)/my
    dz = (z_upper-z_lower)/mz
    dt = 0.90/(material.co*np.sqrt(1.0/(dx**2)+1.0/(dy**2)+1.0/(dz**2)))
    tf = (x_upper-x_lower+5.0)/1.0
    
    return dx,dy,dz,dt,tf

def em2D(mx=128,my=128,mz=128,num_frames=10,cfl=1.0,outdir='./_output',before_step=False,debug=False):
    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    if MPI.COMM_WORLD.rank==0:
        material.dump()
        source.dump()

    num_eqn   = 6
    num_waves = 4
    num_aux   = 12

#   grid pre calculations and domain setup
    dx,dy,dz,dt,tf = grid_basic(x_lower,x_upper,y_lower,y_upper,z_lower,z_upper,mx,my,mz,cfl)
    print tf
    x = pyclaw.Dimension('x',x_lower,x_upper,mx)
    y = pyclaw.Dimension('y',y_lower,y_upper,my)
    z = pyclaw.Dimension('z',z_lower,z_upper,mz)
    
    domain = pyclaw.Domain([x,y,z])

#   Solver settings
    solver = pyclaw.SharpClawSolver3D()
    solver.num_waves  = num_waves
    solver.num_eqn    = num_eqn
    solver.weno_order = 5

    solver.dt_variable = True
    solver.dt_initial  = dt/2.0
    solver.dt_max      = dt
    solver.max_steps   = int(2*tf/dt)
    
#   Import Riemann and Tfluct solvers
    import maxwell_3d_rp
    import maxwell_3d_tfluct

    solver.tfluct_solver = True
    solver.fwave = True
    
    solver.rp = maxwell_3d_rp
    solver.tfluct = maxwell_3d_tfluct
    
    solver.cfl_max = cfl+0.5
    solver.cfl_desired = cfl

#   boundary conditions
    solver.bc_lower[0] = pyclaw.BC.custom
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_lower[2] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap
    solver.bc_upper[2] = pyclaw.BC.extrap
    solver.user_bc_lower = source.scattering_bc

    solver.aux_bc_lower[0]= pyclaw.BC.custom
    solver.aux_bc_lower[1]= pyclaw.BC.extrap
    solver.aux_bc_lower[2]= pyclaw.BC.extrap
    solver.aux_bc_upper[0]= pyclaw.BC.extrap
    solver.aux_bc_upper[1]= pyclaw.BC.extrap
    solver.aux_bc_upper[2]= pyclaw.BC.extrap
    
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
    state.problem_data['co'] = material.co
    state.problem_data['zo'] = material.zo
    state.problem_data['eo'] = material.eo
    state.problem_data['mo'] = material.mo    
    state.problem_data['dx'] = state.grid.x.delta
    state.problem_data['dy'] = state.grid.y.delta
    state.problem_data['dz'] = state.grid.z.delta

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
