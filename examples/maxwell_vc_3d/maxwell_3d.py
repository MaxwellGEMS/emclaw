import numpy as np
from emclaw.utils.materials import Material3D
from emclaw.utils.sources import Source3D
from emclaw.utils import basics

x_lower = 0.0
x_upper = 10.0

y_lower = 0.0
y_upper = 5.0

z_lower = 0.0
z_upper = 5.0

sy = y_upper-y_lower
sx = x_upper-x_lower
sz = z_upper-z_lower

mid_point_y = (y_upper-y_lower)/2.0
mid_point_z = (z_upper-z_lower)/2.0

material = Material3D(shape='homogeneous',metal=False)
material.setup()
material._calculate_n()

source = Source3D(material,shape='off',wavelength=2.0)
source.offset[0] = 5.0
source.transversal_shape = 'plane'
source.heading = 'x'
source.setup()
source.transversal_offset = [sy/2.0,sz/2.0]
source.transversal_width  = [sy/2.0,sz/2.0]

def em3D(mx=64, my=64, mz=64, num_frames=10, cfl=1.0, outdir='./_output', before_step=False, debug=False, nl=False, psi=True):
    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    if np.logical_and(debug, MPI.COMM_WORLD.rank==0):
        material.dump()
        source.dump()

    num_eqn   = 6
    num_waves = 4
    num_aux   = 12

#   grid pre calculations and domain setup
    _, _, _, dt,tf = basics.grid_basic([[x_lower,x_upper,mx], [y_lower,y_upper,my], [z_lower,z_upper,mz]], cfl = cfl, co = material.co, v = source.v)

    x = pyclaw.Dimension(x_lower,x_upper,mx, name = 'x',)
    y = pyclaw.Dimension(y_lower,y_upper,my, name = 'y',)
    z = pyclaw.Dimension(z_lower,z_upper,mz, name = 'z',)
    
    domain = pyclaw.Domain([x, y, z])

#   Solver settings
    solver = pyclaw.SharpClawSolver3D()
    solver.num_waves  = num_waves
    solver.num_eqn    = num_eqn
    solver.lim_type = 2

    solver.dt_variable = True
    solver.dt_initial  = dt/2.0
    solver.dt_max      = dt
    solver.max_steps   = int(2*tf/dt)
    
#   Import Riemann and Tfluct solvers
    from emclaw.riemann import maxwell_3d_nc_rp as maxwell_3d_rp
    from emclaw.riemann import maxwell_3d_nc_tfluct as maxwell_3d_tfluct

    solver.tfluct_solver = True
    solver.fwave = True
    
    solver.rp = maxwell_3d_rp
    solver.tfluct = maxwell_3d_tfluct
    
    solver.cfl_max = cfl+0.5
    solver.cfl_desired = cfl
    solver.reflect_index = [1,0,5]

#   boundary conditions
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_lower[2] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    solver.bc_upper[2] = pyclaw.BC.wall

    solver.aux_bc_lower[0]= pyclaw.BC.wall
    solver.aux_bc_lower[1]= pyclaw.BC.wall
    solver.aux_bc_lower[2]= pyclaw.BC.wall
    solver.aux_bc_upper[0]= pyclaw.BC.wall
    solver.aux_bc_upper[1]= pyclaw.BC.wall
    solver.aux_bc_upper[2]= pyclaw.BC.wall

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
    state.problem_data['nl']     = nl
    state.problem_data['psi']    = psi

    source._delta[0] = state.grid.x.delta
    source._delta[1] = state.grid.y.delta
    source._delta[2] = state.grid.z.delta

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
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em3D)

    import os
    os.system("python ../../emclaw/utils/visualization.py /home/$USER/emclaw/examples/maxwell_vc_3d/_output")
