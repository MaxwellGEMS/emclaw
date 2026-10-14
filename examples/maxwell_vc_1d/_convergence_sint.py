# Import necessary libraries
import sys
import os
import numpy as np

# Adding paths to custom modules
sys.path.append(os.path.realpath('../utils'))
sys.path.append(os.path.realpath('../'))

# Import custom modules
from utils.materials import Material1D
from utils.sources import Source1D

# Defining some constants and initialize material properties
x_lower = 0.
x_upper = 100.

material = Material1D()
material.shape = 'vibrate'
material.bkg_e = 2.0
material.bkg_h = 2.0
material.setup()
material.delta_length = x_upper
material.delta_corner = x_lower
material.delta_omega = 12.0 * np.pi / 200.0

# Create a source object
source = Source1D(material, shape='off', wavelength=2.0)
source.setup()
source.offset = 10.0
source.pulse_width = 2.0

# Define a function for grid calculations
def grid_basic(x_lower, x_upper, mx, cfl):
    dx = (x_upper - x_lower) / mx
    dt = 0.9 * cfl / (material.co * np.sqrt(1.0 / (dx ** 2)))
    tf = 100.0

    return dx, dt, tf

# Define a function for source terms
def dq_source(solver, state, dt):
    """
    Source terms for Maxwell's equations in one spatial dimension.
    Assume aux values are averages on the cell i
    """

    dq = -dt * state.q * state.aux[2:4, :] / state.aux[0:2, :]

    return dq

# Define the main simulation function
def em1D(mx=1024, num_frames=5, cfl=1.0, outdir='./_output', before_step=True, debug=False, chi3=0.0, chi2=0.0, nl=False, psi=True, em=True, homogeneous=True):

    import clawpack.petclaw as pyclaw
    import petsc4py.PETSc as MPI

    # Set material properties if not homogeneous
    if not homogeneous:
        if nl:
            material.chi3_e = chi3
            material.chi2_e = chi2
            if em:
                material.chi3_m = chi3
                material.chi2_m = chi2

    if MPI.COMM_WORLD.rank == 0:
        material._outdir = outdir
        source._outdir = outdir
        material._dump_to_latex()
        source._dump_to_latex()

    num_eqn = 2
    num_waves = 2
    num_aux = 4

    # Grid pre-calculations and domain setup
    dx, dt, tf = grid_basic(x_lower, x_upper, mx, cfl)
    x = pyclaw.Dimension('x', x_lower, x_upper, mx)
    domain = pyclaw.Domain([x])

    # Solver settings
    solver = pyclaw.SharpClawSolver1D()
    solver.num_waves = num_waves
    solver.num_eqn = num_eqn
    solver.weno_order = 5

    solver.dt_variable = True
    solver.dt_initial = dt / 2.0
    solver.dt_max = dt
    solver.max_steps = int(2 * tf / dt)

    # Set the source
    if not psi:
        # Python 2-style print statement
        # print 'using dq_src'
        print('using dq_src')  
        solver.dq_src = dq_source

    # Import Riemann and Tfluct solvers
    if homogeneous:
        import maxwell_1d_rp
    else:
        import maxwell_1d_nl_rp as maxwell_1d_rp

    solver.tfluct_solver = False
    solver.fwave = True

    solver.rp = maxwell_1d_rp

    if solver.tfluct_solver:
        if homogeneous:
            import maxwell_1d_tfluct
        else:
            import maxwell_1d_nl_tfluct as maxwell_1d_tfluct

        solver.tfluct = maxwell_1d_tfluct

    solver.cfl_max = cfl + 0.05
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

    source._dx = state.grid.x.delta
    material._dx = state.grid.x.delta

    # Array initialization
    source.init(state)
    material.init(state)

    state.q = state.q * state.aux[0:2, :]

    # Controller
    claw = pyclaw.Controller()
    claw.tfinal = tf
    claw.num_output_times = num_frames
    claw.solver = solver
    claw.solution = pyclaw.Solution(state, domain)
    claw.outdir = outdir
    claw.write_aux_always = True

    return claw

if __name__ == "__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em1D)