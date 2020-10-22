import sys
import os
import numpy as np
from emclaw.utils.materials import Material2D
from emclaw.utils.sources import Source2D
from emclaw.utils import basics

x_lower = 0.0
x_upper = 50.0

y_lower = 0.0
y_upper = 50.0

sy = y_upper-y_lower
sx = x_upper-x_lower

material = Material2D(shape='homogeneous',metal=False)
material.setup()
material._calculate_n()

def em2D(mx=128,my=128,num_frames=10,cfl=1.0,outdir='./_output',use_petsc=True, before_step=False,debug=False,heading='x',shape='off',nl=False,psi=True,conservative=True):
    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    source = Source2D(material,shape=shape,wavelength=2.0)

    if shape=='off':
        source.offset.fill(5.0)
        if heading=='xy':
            source.offset[0] = sy/2.0
            source.offset[1] = sx/2.0
            tf = 22.0
    else:
        source.offset[0] = -5.0
        source.offset[1] = sy/2.0
        source.transversal_offset = sy/2.0
        source.transversal_width = sy
        source.transversal_shape = 'plane'
    source.setup()
    source.heading = heading
    source.averaged = True

    #   grid pre calculations and domain setup
    _, _, dt,tf = basics.grid_basic([[x_lower,x_upper,mx], [y_lower,y_upper,my]], 
                                    cfl = cfl, co = material.co, v = source.v)

    if (debug and MPI.COMM_WORLD.rank==0):
        material.dump()
        source.dump()

    num_eqn   = 3
    num_waves = 2
    num_aux   = 6

    x = pyclaw.Dimension(x_lower,x_upper,mx, name = 'x')
    y = pyclaw.Dimension(y_lower,y_upper,my, name = 'y')

    domain = pyclaw.Domain([x,y])

#   Solver settings
    solver = pyclaw.SharpClawSolver2D()
    solver.num_waves  = num_waves
    solver.num_eqn    = num_eqn
    solver.reconstruction_order = 5
    solver.lim_type = 2

    solver.dt_variable = True
    solver.dt_initial  = dt/2.0
    solver.dt_max      = dt
    solver.max_steps   = int(2*tf/dt)

#   Import Riemann and Tfluct solvers
    if conservative:
        from emclaw.riemann import maxwell_2d_rp
    else:
        from emclaw.riemann import maxwell_2d_nc_rp as maxwell_2d_rp

    solver.tfluct_solver = True
    solver.fwave = True

    solver.rp = maxwell_2d_rp

    if solver.tfluct_solver:
        if conservative:
            from emclaw.riemann import maxwell_2d_tfluct
        else:
            from emclaw.riemann import maxwell_2d_nc_tfluct as maxwell_2d_tfluct

    solver.tfluct = maxwell_2d_tfluct

    solver.cfl_max = cfl+0.5
    solver.cfl_desired = cfl
    solver.reflect_index = [1,0]

#   boundary conditions
    if shape=='off':
        solver.bc_lower[0] = pyclaw.BC.wall
        solver.aux_bc_lower[0]= pyclaw.BC.wall
    else:
        solver.bc_lower[0] = pyclaw.BC.custom
        solver.aux_bc_lower[0]= pyclaw.BC.custom
        solver.user_bc_lower = source.scattering_bc
        solver.user_aux_bc_lower = material.setaux_lower

    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall

    solver.aux_bc_lower[1]= pyclaw.BC.wall
    solver.aux_bc_upper[0]= pyclaw.BC.wall
    solver.aux_bc_upper[1]= pyclaw.BC.wall

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
    state.problem_data['nl']     = nl
    state.problem_data['psi']    = psi

    source._dx = state.grid.x.delta
    source._dy = state.grid.y.delta

#   array initialization
    source.init(state)
    material.init(state)

    if conservative:
        state.q = state.q*state.aux[0:3,:,:]

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
