#!/usr/bin/env python
# encoding: utf-8

import numpy as np


# -------- GLOBAL SCALAR DEFINITIONS -----------------------------
# ======== all definitions are in m,s,g unit system.
n_frames = 30
# ....... dimensions .............................................
x_lower = 0.0
x_upper = 10.0e-6                    # lenght [m]
y_lower = 0.0
y_upper = 10.0e-6                   # notice that for multilayer this is value will be over-written
z_lower = 0.0
z_upper = 10.0e-6
# ........ material properties ...................................

# vacuum
eo = 8.854187817e-12            # vacuum permittivity   - [F/m]
mo = 4e-7*np.pi                 # vacuum peremeability  - [V.s/A.m]
co = 1.0/np.sqrt(eo*mo)           # vacuum speed of light - [m/s]
zo = np.sqrt(mo/eo)

# material
mat_shape = 'homogenous'           # material definition: homogeneous, interface, rip (moving perturbation), multilayered

# background refractive index and etas
eta = vac = np.ones([6])
bkg_n     = np.ones([3])

bkg_n[0] = np.sqrt(eta[0]*eta[3])
bkg_n[1] = np.sqrt(eta[1]*eta[4])
bkg_n[2] = np.sqrt(eta[2]*eta[5])

vac[0:2] = mo
vac[3:5] = eo

# if interface declare position
x_change = (x_upper - x_lower)/2.0
y_change = (y_upper - y_lower)/2.0
z_change = (z_upper - z_lower)/2.0

# set moving refractive index or gaussian3D parameters
rip_velocity = np.zeros([3,3])
rip_offset   = np.zeros([2,3])
rip_sigma    = np.zeros([2,3])
delta_eta    = np.zeros([3])

rip_offset[0,:].fill((x_upper-x_lower)/2.0)
rip_offset[1,:].fill((y_upper-y_lower)/2.0)
rip_sigma[0,:].fill((x_upper-x_lower)/25.0)
rip_sigma[1,:].fill((y_upper-y_lower)/25.0)
rip_sigma.fill(10e-6)
rip_sigma.fill(10e-6)
rip_sigma2 = rip_sigma**2

delta_eta = np.zeros([3])
delta_eta = 0.1*eta

# set multilayer parameters

# multilayered definition
n_layers = 2
layers = np.zeros([n_layers,9]) # _layer:  eps mu N t chi2e chi2m chi3e chi3m
layers[0,0] = 1.5
layers[0,1] = 1.5
layers[0,2] = 1.5
layers[0,3] = 10
layers[0,4] = 15e-9
layers[1,0] = 2.5
layers[1,1] = 2.5
layers[1,2] = 2.5
layers[1,3] = layers[0,3] - 1
layers[1,4] = 50e-9
N_layers = 5
if mat_shape=='multilayer':
    y_upper = N_layers*np.sum(layers[:,4])+layers[0,4]
    tlp = np.sum(layers[:,4])
    mlp = np.floor(tlp/1e-9)

# set non-linear parameters of the material
chi2 = chi3 = np.zeros( [6] )

# ........ excitation - initial conditoons .......................

# pre-allocate arrays
ex_sigma     = np.ones( [6])    # x,y,t
ex_offset    = np.zeros([6])
ex_amplitude = np.ones( [6])
ex_kvector   = np.zeros([3])

# fill arrays and set respective values
ex_type   = 'off'
ex_lambda = 1e-6
ex_sigma[:] = 1.0*ex_lambda
ex_offset[:]  = (y_upper-y_lower)/2.0

# post calculations
omega    = 2.0*np.pi*co/ex_lambda   # frequency
k        = 2.0*np.pi/ex_lambda      # k vector magnitude
ex_kvector[0:2] = k                   # propagation along the x-direction

# ........ pre-calculations for wave propagation .................

v = co/bkg_n.min()

# Grid - mesh settings
mx = mz = np.floor(10*(x_upper-x_lower)/ex_lambda)
if mat_shape=='multilayer':
    my = np.floor((y_upper-y_lower)/1e-9)
else:
    my = np.floor(10*(y_upper-y_lower)/ex_lambda)

jmx = complex(0,mx)
jmy = complex(0,my)
jmz = complex(0,mz)

ddx = (x_upper-x_lower)/mx
ddy = (y_upper-y_lower)/my
ddz = (z_upper-z_lower)/mz
ddt = dt=0.50/(co*np.sqrt(1.0/(ddx**2)+1.0/(ddy**2)+1.0/(ddz**2)))
max_steps = 1000000
t_final = (x_upper-x_lower)/v 
print t_final

# -------- GLOBAL FUNCTION DEFINITIONS --------------

# refractive index map definition function 
def eta(t,X,Y,Z):
    """
    eta = eta(t,x,y,z)

    This function returns the refractive index map based on general definitions set earlier,
    Gaussian cases support moving RIPs.
    
    x are the coordinate of the grid centers state.grid.e_j.centers, e_j = x 
         aux holds:
         0: epsilon
         1: mu
         2: epsilon_t
         3: mu_t
    """
    z,y,x = np.mgrid[Z.min():Z.max():jmz,Y.min():Y.max():jmy,X.min():X.max():jmx]
    eta_out = np.zeros( [12,len(x),len(y),len(z)], order='F')

    if mat_shape=='gaussian1dx':

        u_d_eta = np.zeros( [3,12,len(x),len(y),len(z)], order='F')
        u_eta_t = np.zeros( [  12,len(x),len(y),len(z)], order='F')
        u_eta   = np.zeros( [  12,len(x),len(y),len(z)], order='F')
        
        for m in range(0,num_dim):
            if m==0:
                d_grid = x
            elif m==1:
                d_grid = y
            elif m==2:
                d_grid = z
            for i in range(0,num_aux):
                u_d_eta[m,i] = d_grid - rip_velocity[m,i]*t - rip_offset[m,i]

        for i in range(0,num_aux):
            u_eta[i]   = sum((u_d_eta[m,i]/rip_sigma[m,i])**2 for m in range(0,num_dim)) 
            u_eta_t[i] = 2.0*sum((rip_velocity[m,i]*u_d_eta[m,i])/(rip_sigma[m,i]**2) for m in range(0,num_dim))
            
            eta_out[i]   = delta_eta[i]*np.exp(-u_eta[i]) + eta[i]
            eta_out[i+6] = u_eta_t[i]*delta_eta[i]*np.exp(-u_eta[i])

    elif mat_shape=='homogeneous':
        for i in range(0,num_aux/2):
            eta_out[i] = eta[i]

    # elif mat_shape=='interfacex':
    #     eta_out[0:num_aux-1,:,:,:] = 1*(x<x_change) + 4*(x>=x_change)

    # elif mat_shape=='interfacey':
    #     eta_out[0:num_aux-1,:,:,:] = 1*(y<y_change/2) + 4*(x>=y_change/2)

    # elif mat_shape=='multilayer':
    #     for n in range(0,N_layers):
    #         yi = n*tlp
    #         for m in range(0,n_layers):
    #             if m==0:
    #                 eta_out[0,:,:,:] = layers[m,0]*(yi<y)*(y<=yi+layers[m,3])
    #                 eta_out[1,:,:,:] = layers[m,1]*(yi<y)*(y<=yi+layers[m,3])
    #                 eta_out[2,:,:,:] = layers[m,2]*(yi<y)*(y<=yi+layers[m,3])
    #             else:
    #                 eta_out[0,:,:,:] = layers[m,0]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                 eta_out[1,:,:,:] = layers[m,1]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])
    #                 eta_out[2,:,:,:] = layers[m,2]*(yi+layers[m-1,3]<y)*(y<=yi+layers[m,3])


    #     eta_out[0,:,:,:] = layers[0,0]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])
    #     eta_out[1,:,:,:] = layers[0,1]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])
    #     eta_out[2,:,:,:] = layers[0,2]*(N_layers*tlp<y)*(y<=N_layers*tlp+layers[0,4])  

    return eta_out

def update_aux(solver,state):
    state.aux = set_aux(state)
    return state

#   next function might be redundant since it already exists as deltan  
def set_aux(state):
    grid = state.grid
    x = grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    y = grid.y.centers_with_ghost(num_ghost)[:num_ghost]
    z = grid.z.centers_with_ghost(num_ghost)[:num_ghost]
    t = state.t
    aux = eta(t,x,y,z)
    return aux

def setaux_lower(state,dim,t,auxbc,num_ghost):
    grid = state.grid
    x = grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    y = grid.y.centers_with_ghost(num_ghost)[:num_ghost]
    z = grid.z.centers_with_ghost(num_ghost)[:num_ghost]
    t = state.t
    
    auxbc[:,:num_ghost,:num_ghost] = eta(t,x,y,z)
    return auxbc

def setaux_upper(state,dim,t,auxbc,num_ghost):
    grid = state.grid
    X = grid.x.centers_with_ghost(num_ghost)[-num_ghost:]
    Y = grid.y.centers_with_ghost(num_ghost)[-num_ghost:]
    Z = grid.z.centers_with_ghost(num_ghost)[-num_ghost:]
    t = state.t
    z,y,x = np.mgrid[z.min():z.max():jmz,y.min():y.max():jmy,x.min():x.max():jmx]
    auxbc[:,-num_ghost:,-num_ghost:,-num_ghost:] = eta(t,x,y,z)
    return auxbc

def scattering_bc(state,dim,t,qbc,num_ghost):
    """
    EM scattering boundary conditions with three components in  TM-mode Ex, Ey, Hz.
    """
    grid = state.grid
    X = grid.x.centers_with_ghost(num_ghost)[:num_ghost]
    Y = grid.y.centers_with_ghost(num_ghost)[:num_ghost]
    Z = grid.Z.centers_with_ghost(num_ghost)[:num_ghost]
    ts = state.t
    z,y,x = np.mgrid[Z.min():Z.max():jmz,Y.min():Y.max():jmy,X.min():X.max():jmx]
    t0 = 0.0
    aux_left_bc = eta(t,X,Y,Z)
    pulseshape  = np.zeros( [len(X),len(Y),len(Z)], order='F')
    harmonic    = np.zeros( [len(X),len(Y),len(Z)], order='F')

    if ex_type=='plane':
        pulseshape = 1.0
        harmonic   = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='gauss-beam':
        pulseshape = np.exp(-(y - ex_yoff)**2/ex_y_sig**2)
        harmonic   = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='gauss_pulse':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2 - (y - ex_yoff - ex_vy*(ts-t0))**2/ex_y_sig**2)
        harmonic   = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='plane_pulse':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2)
        harmonic   = np.sin(ex_kx*x + ex_ky*y - omega*ts)
    elif ex_type=='simple_pulse2D':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2 - (y - ex_yoff - ex_vy*(ts-t0))**2/ex_y_sig**2)
        harmonic   = 1.0
    elif ex_type=='simple_pulse2D_x':
        pulseshape = np.exp(-(x - ex_xoff - ex_vx*(ts-t0))**2/ex_x_sig**2)
        harmonic   = 1.0
    elif ex_type=='off':
        pulseshape = 0.0
        harmonic   = 0.0

    qbc.fill(0.0)
    qbc[2,:num_ghost,:num_ghost,:num_ghost] = zo*ex_amplitude[1]*pulseshape*harmonic
    qbc[4,:num_ghost,:num_ghost,:num_ghost] = ex_amplitude[1]*pulseshape*harmonic

    return qbc

def qinit(state):
    """
    Initial conditions in simulation grid for electromagnetic components q
    """
    
    if ex_type=='off':
        grid = state.grid
        X = grid.x.centers
        Y = grid.y.centers
        Z = grid.z.centers
        z,y,x = np.mgrid[Z.min():Z.max():jmz,Y.min():Y.max():jmy,X.min():X.max():jmx]
        dd1 = (x_upper-x_lower)/5.0
        dd2 = y_upper-y_lower
        sdd = 1e-6
        r2 = (x-dd1/2.0)**2 #+ (y-dd2/2.0)**2
        state.q.fill(0.0)
        state.q[2,:,:,:] = zo*np.exp(-r2/(sdd**2))
        state.q[4,:,:,:] = 1.0*np.exp(-r2/(sdd**2))
    else:
        state.q.fill(0.0)
    
    return state

# -------- MAIN SCRIPT --------------

def em3D(kernel_language='Fortran',before_step=False,iplot=False,htmlplot=False,use_petsc=True,save_outdir='./_test_le1',solver_type='sharpclaw'):

    if use_petsc:
        import clawpack.petclaw as pyclaw
        from petsc4py import PETSc as MPI
    else:
        from clawpack import pyclaw


#   Solver settings

    solver = pyclaw.SharpClawSolver3D()
    solver.num_waves = 6

    solver.lim_type = 2

    solver.dt_initial= ddt
    solver.max_steps = max_steps
    solver.dt_variable = True

    import maxwell_3d
    solver.rp = maxwell_3d

    solver.fwave = True

    solver.cfl_max = 1.5
    solver.cfl_desired = 0.4


#   print some debug information

    # if use_petsc:
    #     if MPI.COMM_WORLD.rank==0:
    #     print 'setup information:'
    #     print 'v_wave=',v
    #     print 'x_lim=',x_upper,' t_f=',t_final 
    #     print 'mx=',mx,'dx=',ddx,'my=',mx,'dy=',ddy, 'dt=',ddt,'N_max=',max_steps
    #     print 'lambda=',ex_lambda,'freq=',omega

    # if before_step:
    #     if MPI.COMM_WORLD.rank==0:
    #         print 'update aux'

        # solver.call_before_step_each_stage = 1
        # solver.before_step = update_aux
    


#   define number of waves (eqn) and aux (eps,mu)
    num_eqn = 6
    num_aux = 12

#   abstract domain and state setup
    x = pyclaw.Dimension('x',x_lower,x_upper,mx)
    y = pyclaw.Dimension('y',y_lower,y_upper,my)
    z = pyclaw.Dimension('z',z_lower,z_upper,mz)

    domain  = pyclaw.Domain([x,y,z])
    
    state   = pyclaw.State(domain,num_eqn,num_aux)
    state.aux = eta(state.t,state.grid.x.centers,state.grid.y.centers,state.grid.z.centers)
    #state.aux = eta(state)


    state.problem_data['dx']    = x.delta
    state.problem_data['dy']    = y.delta
    state.problem_data['dz']    = z.delta
    state.problem_data['chi2']  = chi2
    state.problem_data['chi3']  = chi3
    state.problem_data['eo']    = eo
    state.problem_data['mo']    = mo
    state.problem_data['co']    = co
    state.problem_data['zo']    = zo
    # state.problem_data['vac1']  = mo
    # state.problem_data['vac2']  = eo


#   Boundary conditions
#   solver.user_bc_lower = scattering_bc
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap
    solver.bc_lower[2] = pyclaw.BC.extrap
    solver.bc_upper[2] = pyclaw.BC.extrap

    solver.user_aux_bc_lower = setaux_lower
    solver.user_aux_bc_upper = setaux_upper
    solver.aux_bc_lower[0] = pyclaw.BC.wall
    solver.aux_bc_upper[0] = pyclaw.BC.wall
    solver.aux_bc_lower[1] = pyclaw.BC.wall
    solver.aux_bc_upper[1] = pyclaw.BC.wall
    solver.aux_bc_lower[2] = pyclaw.BC.wall
    solver.aux_bc_upper[2] = pyclaw.BC.wall
#   Initial solution
    qinit(state)


#   controller
    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.tfinal = t_final
    claw.num_output_times = n_frames
    claw.solver = solver
    claw.solution = pyclaw.Solution(state,domain)
    claw.outdir = save_outdir
    claw.write_aux_always = True
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=save_outdir,file_format=claw.output_format)
    if iplot:     pyclaw.plot.interactive_plot(outdir=save_outdir,file_format=claw.output_format)

    return claw


if __name__=="__main__":
    import sys
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(em3D)
    


