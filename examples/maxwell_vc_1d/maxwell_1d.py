# Import necessary modules
import numpy as np
from emclaw.utils.materials import Material1D
from emclaw.utils.sources import Source1D
from emclaw.utils import basics

# Create a 1D material with a homogeneous shape and set it up
material = Material1D(shape='homogeneous')
material.setup()

# Create a 1D source with shape 'off', wavelength 2.0, and set it up
source = Source1D(material, shape='off', wavelength=2.0)
source.offset = 25.0
source.setup()

x_lower = 0.
x_upper = 100.

# Define a function em1D with various parameters
def em1D(mx=1024, num_frames=10, use_petsc=True, cfl=1.0, conservative=True,
         chi3=0.0, chi2=0.0, nl=False, psi=True, em=True, before_step=False,
         debug=False, outdir='./_output', output_style=1):

    # Check if petsc4py is used, and import necessary modules
    if use_petsc:
        import clawpack.petclaw as pyclaw
        import petsc4py.PETSc as MPI
    else:
        from clawpack import pyclaw

    if nl:
        material.chi3_e = chi3
        material.chi2_e = chi2
        if em:
            material.chi3_m = chi3
            material.chi2_m = chi2

    # Set output directories based on MPI rank if using petsc4py
    if np.logical_and(use_petsc, MPI.COMM_WORLD.rank == 0):
        basics.set_outdirs(material, source, outdir=outdir, debug=debug)
    else:
        basics.set_outdirs(material, source, outdir=outdir, debug=debug)

    num_eqn = 2
    num_waves = 2
    num_aux = 4

    # Grid precalculations and domain setup
    dx, dt, tf = basics.grid_basic([[x_lower, x_upper, mx]], cfl, material.co, source.v)
    x = pyclaw.Dimension(x_lower, x_upper, mx, name='x')
    domain = pyclaw.Domain([x])

    # Solver settings
    solver = pyclaw.SharpClawSolver1D()
    solver.num_waves = num_waves
    solver.num_eqn = num_eqn
    solver.lim_type = 2

    solver.dt_variable = True
    solver.dt_initial = dt / 2.0
    solver.dt_max = dt
    solver.max_steps = int(2 * tf / dt)
    
    # Import Riemann and Tfluct solvers based on 'conservative' flag
    if conservative:
        from emclaw.riemann import maxwell_1d_rp
    else:
        from emclaw.riemann import maxwell_1d_nc_rp as maxwell_1d_rp

    solver.tfluct_solver = True
    solver.fwave = True
    
    solver.rp = maxwell_1d_rp

    if solver.tfluct_solver:
        if conservative:
            from emclaw.riemann import maxwell_1d_tfluct
        else:
            from emclaw.riemann import maxwell_1d_nc_tfluct as maxwell_1d_tfluct

        solver.tfluct = maxwell_1d_tfluct
    
    solver.cfl_max = cfl + 0.5
    solver.cfl_desired = cfl

    # Boundary conditions
    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall

    solver.aux_bc_lower[0] = pyclaw.BC.wall
    solver.aux_bc_upper[0] = pyclaw.BC.wall

    solver.reflect_index = [0]

    # Before step configure
    if before_step:
        solver.call_before_step_each_stage = True
        solver.before_step = material.update_aux

    # State setup
    state = pyclaw.State(domain, num_eqn, num_aux)
    
    state.problem_data['chi2_e'] = material.chi2_e
    state.problem_data['chi3_e'] = material.chi3_e
    state.problem_data['chi2_m'] = material.chi2_m
    state.problem_data['chi3_m'] = material.chi3_m
    state.problem_data['eo'] = material.eo
    state.problem_data['mo'] = material.mo
    state.problem_data['co'] = material.co
    state.problem_data['zo'] = material.zo
    state.problem_data['dx'] = state.grid.x.delta
    state.problem_data['nl'] = nl
    state.problem_data['psi'] = psi
    state.problem_data['conservative'] = conservative

    source._dx = state.grid.x.delta
    material._dx = state.grid.x.delta

    # Array initialization
    source.init(state)
    material.init(state)

    # Multiply state.q by state.aux[0:2,:] if 'conservative' is True
    if conservative:
        state.q = state.q * state.aux[0:2, :]

    # Controller
    claw = pyclaw.Controller()
    claw.tfinal = tf
    claw.num_output_times = num_frames
    claw.solver = solver
    claw.solution = pyclaw.Solution(state, domain)
    claw.outdir = outdir
    claw.write_aux_always = True
    claw.output_style = output_style
    claw.num_output_times = 10

    return claw

if __name__ == "__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em1D)
