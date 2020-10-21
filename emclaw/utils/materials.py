import numpy as np
import pickle
import os

dim = {}
dim['xyz'] = [0,1,2]
dim['xy' ] = [0,1]
dim['xz' ] = [0,2]
dim['yz' ] = [1,2]
dim['x'  ] = [0]
dim['y'  ] = [1]
dim['z'  ] = [2]

def user_material():
    pass

class Material(object):

    def save(self):
        print(self._outdir)
        """save class as self.name.txt"""
        try:
            os.makedirs(self._outdir)
        except:
            pass
        d = self
        d.function = None
        d.delta_function = None
        d.delta_smooth_function = None
        f = open(os.path.join(self._outdir,self.name + '.saved'),'w')
        pickle.dump(d, f)
        f.close()

    def _set_vacuum(self):
        if self.normalized_vacuum:
            self.eo = 1.0
            self.mo = 1.0
        else:
            self.eo = 8.854187817e-12
            self.mo = 4e-7*np.pi

        self.co = 1.0/np.sqrt(self.eo*self.mo)
        self.zo = np.sqrt(self.mo/self.eo)

        return self.eo,self.mo

    def _unpack_options(self,options={}):
        # first add the options passed to the function
        for key in options:
            setattr(self,key,options[key])

        # unpack options from self.options{}
        for key in self.options:
            setattr(self,key,self.options[key])
    
    def dump(self):
        for attr in sorted(dir(self)):
            if not attr.startswith('_'):
                print("%s = %s" % (attr, getattr(self, attr)))

    def _dump_to_latex(self):
        from tabulate import tabulate
        strt = r'\begin{table}[h!]' + '\n' + r'\centering' + '\n' + r'\begin{tabular}[cl]' + '\n' + r'\hline' + '\n'
        strt = strt + r'variable & value(s) \\' + '\n' + r'\hline' +'\n'
        for attr in sorted(dir(self)):
            if not attr.startswith('_'):
                s = getattr(self,attr)
                if isinstance(s, str):
                    strt = strt + '\t' + r'\verb+' + attr + '+ \t' + r'&' + '\t' + s + r' \\' + '\n'
                elif isinstance(s,float):
                    strt = strt + '\t' + r'\verb+' + attr + '+ \t' + r'&' + '\t' + str(s) + r' \\' + '\n'
                elif isinstance(s,bool):
                    strt = strt + '\t' + r'\verb+' + attr + '+ \t' + r'&' + '\t' + str(s) + r' \\' + '\n'
                else:
                    try:
                        len(s)
                        strt = strt + '\t' + r'\multicolumn{1}{c}\multirow{'+str(np.shape(s)[0])+r'}{*}{\verb+' + attr + r'+}' + '\t' + r'&' + '\t'
                        for k in range(np.shape(s)[0]):
                            strt = strt + str(s[k]) + r' \\'
                        strt = strt + '\n'
                    except:
                        if ('function' in str(s)): s=str(s).split('function ')[1].split('at')[0]
                        if ('method' in str(s)): s=str(s).split('method')[1].split('.')[1].split('of')[0]
                        if ('ufunc' in str(s)): s=str(s).split('ufunc ')[1].split('>')[0]
                        strt = strt + '\t' + r'\verb+' + attr + '+ \t' + r'&' + '\t' + str(s) + r' \\' + '\n'
        strt = strt + r'\end{tabular}' + '\n' + r'\end{table}' + '\n'
        import uuid
        try:
            os.makedirs(self._outdir)
        except:
            pass
        f = open(os.path.join(self._outdir,'_material_'+str(uuid.uuid1())+'.tex'),'a')
        f.write(strt)
        f.close()

    def _dump(self,obj):
        for attr in sorted(dir(obj)):
            try:
                print("%s = %s" % (attr, getattr(obj, attr)))
            except:
                pass

    def plot(self,eta):
        import matplotlib
        # set matplotlib to work over X-forwarding
        matplotlib.use('Agg')
        from matplotlib import pylab as plt
        plt.figure()
        plt.pcolormesh(eta)
        plt.draw()
        plt.savefig('./debug.png',dpi=320)

    def setaux_lower(self,state,dim,t,qbc,auxbc,num_ghost):
        grid = state.grid
        grid.compute_c_centers_with_ghost(num_ghost,recompute=True)
        
        if state.num_dim==1:
            x = grid.x.centers_with_ghost[:num_ghost]
            auxbc[:,:num_ghost] = self.function(x,t)
        elif state.num_dim==2:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            if dim.name==state.grid.dimensions[0].name:
                x = x[:num_ghost,:]
                y = y[:num_ghost,:]
                auxbc[:,:num_ghost,:] = self.function(x,y,t)
            else:
                x = x[:,:num_ghost]
                y = y[:,:num_ghost]
                auxbc[:,:,:num_ghost] = self.function(x,y,t)
        elif state.num_dim==3:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            z = grid._c_centers_with_ghost[2]
            if dim.name==state.grid.dimensions[0].name:
                x = x[:num_ghost,:,:]
                y = y[:num_ghost,:,:]
                z = z[:num_ghost,:,:]
                auxbc[:,:num_ghost,:,:] = self.function(x,y,z,t)
            elif dim.name==state.grid.dimensions[1].name:
                x = x[:,:num_ghost,:]
                y = y[:,:num_ghost,:]
                z = z[:,:num_ghost,:]
                auxbc[:,:,:num_ghost,:] = self.function(x,y,z,t)
            elif dim.name==state.grid.dimensions[2].name:
                x = x[:,:,:num_ghost]
                y = y[:,:,:num_ghost]
                z = z[:,:,:num_ghost]
                auxbc[:,:,:,:num_ghost] = self.function(x,y,z,t)  

        return auxbc

    def setaux_upper(self,state,dim,t,qbc,auxbc,num_ghost):
        grid = state.grid
        grid.compute_c_centers_with_ghost(num_ghost,recompute=True)
        if state.num_dim==1:
            x = grid.x.centers_with_ghost[-num_ghost:]
            auxbc[:,-num_ghost:] = self.function(x,t)
        elif state.num_dim==2:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            if dim.name==state.grid.dimensions[0].name:
                x = x[-num_ghost:,:]
                y = y[-num_ghost:,:]
                auxbc[:,-num_ghost:,:] = self.function(x,y,t)
            else:
                x = x[:,-num_ghost:]
                y = y[:,-num_ghost:]
                auxbc[:,:,-num_ghost:] = self.function(x,y,t)
        elif state.num_dim==3:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            z = grid._c_centers_with_ghost[2]
            if dim.name==state.grid.dimensions[0].name:
                x = x[-num_ghost:,:,:]
                y = y[-num_ghost:,:,:]
                z = z[-num_ghost:,:,:]
                auxbc[:,-num_ghost:,:,:] = self.function(x,y,z,t)
            elif dim.name==state.grid.dimensions[1].name:
                x = x[:,-num_ghost:,:]
                y = y[:,-num_ghost:,:]
                z = z[:,-num_ghost:,:]
                auxbc[:,:,-num_ghost:,:] = self.function(x,y,z,t)
            elif dim.name==state.grid.dimensions[2].name:
                x = x[:,:,-num_ghost:]
                y = y[:,:,-num_ghost:]
                z = z[:,:,-num_ghost:]
                auxbc[:,:,:,-num_ghost:] = self.function(x,y,z,t)  

        return auxbc

    def update_aux(self,solver,state):
        grid = state.grid
        grid.compute_c_centers()
        t = state.t
        if state.num_dim==1:
            x = grid.x.centers
            state.aux = self.function(x,t)

        if state.num_dim==2:
            x,y = grid._c_centers
            state.aux = self.function(x,y,t)

        if state.num_dim==3:
            x,y,z = grid._c_centers
            state.aux = self.function(x,y,z,t)

        return state

    def impose_metal(self,solver,state):
        if self.update_at_each_stage:
            self.update_aux(solver,state)
        if state.num_dim==2:
            grid = state.grid
            x,y = grid.c_centers
            for k in range(0,len(self.metal_corners)):
                xi1,xi2 = self.metal_corners[k,:,0]
                yi1,yi2 = self.metal_corners[k,:,1]
                state.q[0:2,(x>=xi1)*(x<=xi2)*(y>=yi1)*(y<=yi2)] = 0.0
        return state

    def init(self,state):
        grid = state.grid
        grid._compute_c_centers()
        t = state.t
        if state.num_dim==1:
            x = grid.x.centers
            state.aux = self.function(x,t)
        if state.num_dim==2:
            x,y = grid._c_centers
            state.aux = self.function(x,y,t)
        if state.num_dim==3:
            x,y,z = grid._c_centers
            state.aux = self.function(x,y,z,t)
        return state

    def _get_vibrate(self,x=0.0,y=0.0,span=None,t=0.0):

        w  = self.delta_function(self.delta_omega*t)
        dw = self.delta_sign_dt*self.delta_omega*self.delta_function_dt(self.delta_omega*t)

        if self.delta_smooth:
            s  = self.delta_smooth_function(x,y)
            w  = s*w
            dw = s*dw

        if span is not None:
            w  =  w*span
            dw = dw*span

        return w,dw

    def general_setup(self,options={}):
        self._unpack_options(options=options)
        self._set_vacuum()

        temp_flag = False
        if self.custom:
            self.custom_function = user_material

        if self.shape=='homogeneous':
            self.function = self._homogeneous
        
        if self.shape.startswith('moving'):
            self.offset         = np.zeros([self.num_dim,self.num_aux/2])
            self.delta_velocity = np.zeros([self.num_dim,self.num_aux/2])
            self._moving        = True
            self.update_at_each_stage = True

            if 'gauss' in self.shape:
                self.function = self._gaussian_rip
            
            if 'tanh' in self.shape:
                self.function = self._tanh_rip

            self.offset[0,:].fill(10.0)
            self.delta_velocity[0,:].fill(0.59)
        
        if 'gauss' in self.shape:
            temp_flag = True

        if 'tanh' in self.shape:
            temp_flag = True

        self.temp_flag = temp_flag
        if temp_flag:
            self.delta_sigma        = 5.0*np.ones([self.num_dim,self.num_aux/2])
            self.relative_amplitude = 0.1*np.ones([self.num_aux/2])
            self.delta_eta          = self.relative_amplitude*self.bkg_eta
            self.em_equal           = True
            self.delta_sign         = 1.0

            if not self._moving:
                self.function   = self._gaussian
                self.delta_sign = 1.0
            
            self._rip_precalc = False

        if self.shape.startswith('fiber'):
            self.fiber_eta = np.ones([self.num_aux/2])
            self.fiber_corner = np.zeros([self.num_dim])
            self.fiber_corner[0] = -5.0
            self.fiber_width  = 5.0
            self.fiber_length = 100.0

        if self.shape=='fiber single':
            self.function     = self._single_fiber

        if self.shape=='fiber double':
            self.fiber_eta    = np.ones([2,self.num_aux/2])
            self.fiber_corner = np.zeros([2,self.num_dim])
            self.fiber_width  = 5.0*np.ones([self.num_aux/2])
            self.fiber_length = 100.0*np.ones([self.num_aux/2])
            self.function     = self._double_fiber

        if self.shape=='fiber vibrate':
            self.delta_width  = 5.0
            self.delta_length = 5.0
            self.delta_corner = np.zeros([self.num_dim])
            self.delta_corner[0] = 5.0
            self.delta_eta    = np.ones([self.num_aux/2])
            self.delta_smooth = False
            self.delta_omega  = 2.0*np.pi
            self.delta_function    = np.cos
            self.delta_function_dt = np.sin
            self.delta_sign_dt     = -1.0
            self.delta_angular_velocity = None
            self.delta_smooth_function = self._gaussianf
            self.delta_smooth_width  = 5.0
            self.delta_smooth_length = 5.0
            self.delta_smooth_np     = 2.0*np.ones([self.num_dim])
            self.function = self._oscillate
            
            self.update_at_each_stage = True
        
        if self.shape=='expansion':
            self.delta_radii          = 5.0*np.ones([self.num_aux/2])
            self.delta_expansion_rate = np.zeros([self.num_aux/2])
            self.delta_velocity       = np.zeros([self.num_dim,self.num_aux/2])
            self.offset               = np.zeros([self.num_dim,self.num_aux/2])
            self.delta_sigma          = np.ones([self.num_dim,self.num_aux/2])
            self.delta_eta            = 0.1*np.ones([self.num_aux/2])
            self.update_at_each_stage = True
            self.delta_sign           = 1.0

            self.offset[0,:].fill(10.0)
            self.dim = 'xy'
            self.function = self._expanding

        if self.nonlinear:
            self.chi2 = np.zeros( [self.num_aux/2], order='F')
            self.chi3 = np.zeros( [self.num_aux/2], order='F')

        if self.metal:
            self.metal_corners = []

        return

    def __init__(self,normalized=True,shape='homogeneous'):
        self.normalized_vacuum = normalized
        self.shape     = shape
        self.custom    = False
        self.averaged  = False
        self.nonlinear = True
        self.bkg_eta   = np.ones([self.num_aux//2])
        self.delta     = np.ones([self.num_aux//2])
        self.bkg_n     = np.ones([self.num_dim])
        self.n_max     = np.ones([self.num_dim])
        self._outdir   = './_output'
        self.name      = 'material'
        self._moving   = False
        self.v         = 1.0
        self.update_at_each_stage = False
        self.function  = None

class Material1D(Material):

    def setup(self,options={}):
        self._unpack_options(options=options)
        self._set_vacuum()

        self.bkg_n = np.sqrt(self.bkg_e*self.bkg_h)
        self.n_max = self.bkg_n

        temp_flag = False

        if self.custom:
            self.custom_function = user_material

        if self.shape=='homogeneous':
            self.function = self._homogeneous
        
        if self.shape.startswith('moving'):
            self.delta_velocity_e = 0.59
            self.delta_velocity_h = self.delta_velocity_e
            self.offset_e   = 10.0
            self.offset_h   = self.offset_e
            self._moving    = True

            if 'gauss' in self.shape:
                self.function = self._gaussian_rip
            
            if 'tanh' in self.shape:
                self.function = self._tanh_rip
        
        if 'gauss' in self.shape:
            temp_flag = True

        if 'tanh' in self.shape:
            temp_flag = True

        self.temp_flag = temp_flag
        if temp_flag:
            self.sigma_e            = 5.0
            self.sigma_h            = self.sigma_e
            self.relative_amplitude = 0.1
            self.delta_n     = self.relative_amplitude*(self.bkg_n)
            self.delta_e     = self.delta_n
            self.delta_h     = self.delta_n
            self.em_equal = True

            if not self._moving:
                self.function = self._gaussian
            
            self._rip_precalc = False

        if self.shape=='vibrate':
            self.delta_length = 5.0
            self.delta_corner = 5.0
            self.delta_e      = 1.0
            self.delta_h      = 1.0
            self.delta_smooth = False
            self.delta_omega  = 2.0*np.pi
            self.delta_function    = np.cos
            self.delta_function_dt = np.sin
            self.delta_sign_dt     = -1.0
            self.delta_smooth_function = self._gaussianf
            self.function = self._oscillate

        if self.shape=='farago':
            self.function = self._farago

        if self.nonlinear:
            self.chi2_e = 0.0
            self.chi3_e = 0.0
            self.chi2_m = 0.0
            self.chi3_m = 0.0

        return

    def _farago(self,x,t):
        eta = np.zeros( [4,len(x)], order='F')

        eta[0,:].fill((1.0+t)**2)
        eta[1,:].fill(1.0)
        eta[2,:].fill(2.0*(1.0+t))
        eta[3,:].fill(0.0)

        return eta

    def _x_w_offset(self,x,v=[0.0,0.0],t=0.0):
        u_x_e = x - v[0]*t - self.offset_e
        u_x_m = x - v[1]*t - self.offset_h

        return u_x_e,u_x_m

    def _gaussian_rip(self,x,t):
                
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e,u_x_m = self._x_w_offset(x,v=[self.delta_velocity_e,self.delta_velocity_h],t=t)

        u_e_t = 2.0*((self.delta_velocity_e*u_x_e)/(self.sigma_e**2))
        u_m_t = 2.0*((self.delta_velocity_h*u_x_m)/(self.sigma_h**2))
        u_e = (u_x_e/self.sigma_e)**2
        u_m = (u_x_m/self.sigma_h)**2

        if self.averaged:
            from scipy.special import erf
            ddx = self._dx/2.0
            arg1_e = (ddx - u_x_e)/self.sigma_e
            arg2_e = (ddx + u_x_e)/self.sigma_e
            arg1_m = (ddx - u_x_m)/self.sigma_h
            arg2_m = (ddx + u_x_m)/self.sigma_h
            eta[0,:] = (1/self._dx)*(2.0*ddx*self.bkg_e + 0.5*self.delta_e*np.sqrt(np.pi)*self.sigma_e*(erf(arg1_e)+erf(arg2_e)))
            eta[1,:] = (1/self._dx)*(2.0*ddx*self.bkg_h + 0.5*self.delta_h*np.sqrt(np.pi)*self.sigma_h*(erf(arg1_m)+erf(arg2_m)))
        else:
            eta[0,:] = self.delta_e*np.exp(-u_e) + self.bkg_e
            eta[1,:] = self.delta_h*np.exp(-u_m) + self.bkg_h

        eta[2,:] = u_e_t*self.delta_e*np.exp(-u_e)
        eta[3,:] = u_m_t*self.delta_h*np.exp(-u_m)

        return eta

    def _gaussian(self,x,t=0):
        
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e,u_x_m = self._x_w_offset(x,)

        u_e = (u_x_e/selfsigma_e)**2
        u_m = (u_x_m/selfsigma_h)**2
        eta[0,:] = self.delta_e*np.exp(-u_e) + self.bkg_e
        eta[1,:] = self.delta_h*np.exp(-u_m) + self.bkg_h    

        return eta

    def _gaussianf(self,x):

        u = x - self.delta_corner + self.delta_length/2.0

        r2 = u**2/self.delta_length**2
        
        g = np.exp(-r2)
        
        return g

    def _oscillate(self,x,t=0):

        eta = np.zeros( [4,x.shape[0]], order='F')

        xid = self.delta_corner
        
        spand = ((x>=xid)*(x<=(xid+self.delta_length)))

        w,dw = self._get_vibrate(span=spand,t=t)

        eta[0,:] = self.bkg_e + self.delta_e*w
        eta[1,:] = self.bkg_h + self.delta_h*w
        eta[2,:] = self.delta_e*dw
        eta[3,:] = self.delta_h*dw

        return eta

    def _tanh_rip(self,x,t):
        
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e,u_x_m = self._x_w_offset(x,v=[self.delta_velocity_e,self.delta_velocity_h],t=t)

        eta[0,:] = (self.delta_e/2.0)*(1.0 + np.tanh(u_x_e)) + self.bkg_e
        eta[1,:] = (self.delta_h/2.0)*(1.0 + np.tanh(u_x_m)) + self.bkg_h
        eta[2,:] = -(self.delta_e*self.delta_velocity_e/(2.0*self.sigma_e))/(np.cosh(u_x_e)**2)
        eta[3,:] = -(self.delta_h*self.delta_velocity_h/(2.0*self.sigma_h))/(np.cosh(u_x_m)**2)

        return eta

    def _homogeneous(self,x,t=0):

        eta = np.zeros( [4,len(x)], order='F')

        eta[0,:] = self.bkg_e
        eta[1,:] = self.bkg_h

        return eta

    def _calculate_n(self):
        
        eta = self.bkg_eta
        
        if hasattr(self,'fiber_eta'):
            eta = eta + self.fiber_eta

        if hasattr(self,'delta_eta'):
            eta = eta + self.delta_eta

        self.bkg_n = np.sqrt(self.bkg_eta[0]*self.bkg_eta[2])
        self.n_max = np.sqrt(eta[0]*eta[2])

        return

    def __init__(self,normalized=True,shape='homogeneous'):
        self.num_aux = 4
        self.num_dim = 1
        self.options = {}
        super(Material1D,self).__init__(normalized,shape)
        self.dim = 'x'
        self.bkg_e = 1.0
        self.bkg_h = 1.0
        self._dx = 1

class Material2D(Material):

    def setup(self,options={}):
        self.general_setup()

    def _calculate_n(self):
        
        eta = self.bkg_eta
        
        if hasattr(self,'fiber_eta'):
            if len(self.fiber_eta)==1:
                eta = eta + self.fiber_eta
            else:
                eta = eta + self.fiber_eta[0]

        if hasattr(self,'delta_eta'):
            eta = eta + self.delta_eta

        self.bkg_n[0] = np.sqrt(self.bkg_eta[0]*self.bkg_eta[2])
        self.bkg_n[1] = np.sqrt(self.bkg_eta[1]*self.bkg_eta[2])

        self.n_max[0] = np.sqrt(eta[0]*eta[2])
        self.n_max[1] = np.sqrt(eta[1]*eta[2])

        return

    def _gaussian_rip(self,x,y,t):
                
        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')
        _r2 = np.zeros( [3,x.shape[0],y.shape[1]], order='F')
        _rt = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        u = v = False

        if self.dim=='x': u = True
        
        if self.dim=='y': v =  True
        
        if self.dim=='xy': u = v = True

        if u:
            for i in range(0,3):
                _temp1 = x - self.offset[0,i] - self.delta_velocity[0,i]*t
                _r2[i] = (_temp1/self.delta_sigma[0,i])**2
                _rt[i] = (_temp1*self.delta_velocity[0,i])/(self.delta_sigma[0,i]**2)

        if v:
            for i in range(0,3):
                _temp2 = y - self.offset[1,i] - self.delta_velocity[1,i]*t
                _r2[i] = _r2[i] + (_temp2/self.delta_sigma[1,i])**2
                _rt[i] = _rt[i] + (_temp2*self.delta_velocity[1,i])/(self.delta_sigma[1,i]**2)

        _r2 = np.exp(-_r2)
        _rt = 2.0*_r2*_rt

        for i in range(0,3):
            eta[i  ] = self.delta_sign*self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+3] = self.delta_sign*self.delta_eta[i]*_rt[i]

        return eta

    def _averaged_gauss(self,x,dx=None,s=1.0,xo=0.0,v=0.0,t=0.0):
        from scipy.special import erf
        arg = xo + v*t - x
        if dx is None:
            dx = self._dx
        ddx = dx/2.0

        erravg = (np.sqrt(np.pi)*s*(erf((ddx + arg)/s) + erf((ddx - arg)/s)))/(2.0*dx)
        
        return erravg

    def _gaussian_rip_averaged(self,x,y,t=0):
        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')
        _r2 = np.ones(  [3,x.shape[0],y.shape[1]], order='F')
        _rp = np.ones(  [3,x.shape[0],y.shape[1]], order='F')
        _rt = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        u = v = False

        if self.dim=='x': u = True

        if self.dim=='y': v =  True

        if self.dim=='xy': u = v = True

        if u:
            for i in range(0,3):
                _r2[i] = self._averaged_gauss(x,s=self.delta_sigma[0,i],xo=self.offset[0,i],v=self.delta_velocity[0,i],t=t)
                _temp1 = x - self.offset[0,i] - self.delta_velocity[0,i]*t
                _rp[i] = (_temp1/self.delta_sigma[0,i])**2
                _rt[i] = (_temp1*self.delta_velocity[0,i])/(self.delta_sigma[0,i]**2)

        if v:
            for i in range(0,3):
                _r2[i] = (_r2[i])*self._averaged_gauss(y,s=self.delta_sigma[1,i],xo=self.offset[1,i],v=self.delta_velocity[1,i],t=t)
                _temp2 = y - self.offset[1,i] - self.delta_velocity[1,i]*t
                _rp[i] = _rp[i] + (_temp2/self.delta_sigma[1,i])**2
                _rt[i] = _rt[i] + (_temp2*self.delta_velocity[1,i])/(self.delta_sigma[1,i]**2)

        _rp = np.exp(-_rp)
        _rt = 2.0*_rp*_rt

        for i in range(0,3):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+3] = self.delta_eta[i]*_rt[i]

        return eta

    def _gaussianf(self,x,y):

        u = x - (self.delta_corner[0] + self.delta_smooth_length/2.0)
        v = y - (self.delta_corner[1] + self.delta_smooth_width/2.0)

        r2 = u**2/self.delta_length**2 + v**2/self.delta_width**2
        
        g = np.exp(-r2)
        
        return g

    def _sinsin(self,x,y):

        l = self.delta_smooth_length
        w = self.delta_smooth_width

        s = np.sin(self.delta_smooth_np[0]*x*np.pi/l)*np.sin(self.delta_smooth_np[1]*y*np.pi/w)
        s = s*(((x/l)<=1.0)*((y/w)<=1.0))

        return s

    def _rotating_sinsin(self,x,y,t):

        p = self.delta_angular_velocity
        o = self.delta_omega
        l = self.delta_smooth_length
        w = self.delta_smooth_width
        xo = l/2
        yo = w/2

        x = self.delta_smooth_np[0]*(x-xo)*np.pi/xo
        y = self.delta_smooth_np[1]*(y-yo)*np.pi/yo

        s  = np.sin(o*t)*(np.sin(x*np.cos(p*t) - y*np.sin(p*t))*np.sin(x*np.sin(p*t) + y*np.cos(p*t)))

        ds = -p*np.cos(x*np.cos(p*t) - y*np.sin(p*t))*np.sin(o*t)*(y*np.cos(p*t) + x*np.sin(p*t))* \
            np.sin(y*np.cos(p*t) + x*np.sin(p*t)) + (p*np.cos(y*np.cos(p*t) + x*np.sin(p*t))* \
            np.sin(o*t)*(x*np.cos(p*t) - y*np.sin(p*t)) + o*np.cos(o*t)*np.sin(y*np.cos(p*t) + x*np.sin(p*t)))* \
            np.sin(x*np.cos(p*t) - y*np.sin(p*t))
        
        return s,ds

    def _gaussian(self,x,y,t=0):
        
        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')
        _r2 = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        u = v = False

        if self.dim=='x': u = True
        
        if self.dim=='y': v =  True
        
        if self.dim=='xy': u = v = True

        if u:
            for i in range(0,3): _r2[i] = ((x - self.offset[0,i])/self.sigma[0,i])**2

        if v:
            for i in range(0,3): _r2[i] = _r2[i] + ((y - self.offset[0,i])/self.sigma[0,i])**2

        for i in range(0,3): eta[i] = self.delta_eta[i]*np.exp(-_r2[i]) + self.bkg_eta[i]

        return eta

    def _tanh_rip(self,x,t):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')
        _r2 = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        u = v = False
        
        if self.dim=='x': u = True
        
        if self.dim=='y': v = True

        if self.dim=='xy': u = v = True

        if u:
            for i in range(0,3): _r2[i] = (x - self.delta_velocity[0,i] - self.offset[0,i])/self.delta_sigma[0,i]
            _rt = self.delta_velocity[0,:]/(2.0*self.delta_sigma[0,:])
        if v:
            for i in range(0,3): _r2[i] = _r2[i] + (y - self.delta_velocity[1,i] - self.offset[1,i])/self.delta_sigma[1,i]
            _rt = _rt + self.delta_velocity[1,:]/(2.0*self.delta_sigma[1,:])

        for i in range(0,3):
            eta[i]   = (self.delta_eta[i]/2.0)*(1.0 + np.tanh(_r2[i])) + self.bkg_eta[i]
            eta[i+3] = -(self.delta_eta[i]*_rt[i])*(1.0/np.cosh(_r2[i])**2)

        return eta

    def _homogeneous(self,x,y,t=0):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')

        for i in range(0,3): eta[i] = self.bkg_eta[i]

        return eta

    def _single_fiber(self,x,y,t=0):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')

        yi = self.fiber_corner[1]
        xi = self.fiber_corner[0]

        spany = (y>=yi)*(y<=(yi + self.fiber_width ))
        spanx = (x>=xi)*(y<=(xi + self.fiber_length))

        span  = spanx*spany

        for i in range(0,3): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span

        return eta

    def _double_fiber(self,x,y,t=0):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')

        y1i,y2i = self.fiber_corner[:,1]
        x1i,x2i = self.fiber_corner[:,0]

        span1y = ((y>=y1i)*(y<=(y1i+self.fiber_width[0])))
        span2y = ((y>=y2i)*(y<=(y2i+self.fiber_width[1])))

        span1x = ((x>=x1i)*(x<=(x1i+self.fiber_length[0])))
        span2x = ((x>=x2i)*(x<=(x2i+self.fiber_length[1])))

        span1 = span1y*span1x
        span2 = span2y*span2x

        for i in range(0,3): eta[i] = self.bkg_eta[i] + self.fiber_eta[0,i]*span1 + self.fiber_eta[1,i]*span2

        return eta

    def _oscillate(self,x,y,t=0):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')

        xi,yi = self.fiber_corner

        span = ((y>=yi)*(y<=(yi+self.fiber_width)))*((x>=xi)*(x<=(xi+self.fiber_length)))

        xid,yid = self.delta_corner
        
        spand = ((y>=yid)*(y<=(yid+self.delta_width)))*((x>=xid)*(x<=(xid+self.delta_length)))

        if self.delta_angular_velocity is None:
            w,dw = self._get_vibrate(x=x,y=y,span=spand,t=t)
        else:
            w,dw = self._rotating_sinsin(x,y,t)
            w  = w*spand
            dw = dw*spand

        for i in range(0,3): 
            eta[i]   = self.bkg_eta[i] + self.fiber_eta[i]*span + self.delta_eta[i]*w
            eta[i+3] = self.delta_eta[i]*dw

        return eta

    def _expanding(self,x,y,t=0.0):
        for k in range(2): self.delta_sigma[k,:].fill(self.delta_radii[k]/((1+self.delta_expansion_rate[k]*t)))

        eta = self._gaussian_rip(x,y,t)

        return eta

    def __init__(self,normalized=True,shape='homogeneous',metal=False):
        self.num_aux = 6
        self.num_dim = 2
        self.options = {}
        super(Material2D,self).__init__(normalized,shape)
        self.dim = 'x'
        self.metal = metal
        self._dx = 1.0
        self._dy = 1.0

class Material3D(Material):
    def setup(self,options={}):
        self.general_setup()

    def _calculate_n(self):
        
        eta = self.bkg_eta
        
        if hasattr(self,'fiber_eta'):
            eta = eta + self.fiber_eta

        if hasattr(self,'delta_eta'):
            eta = eta + self.delta_eta

        self.bkg_n[0] = np.sqrt(self.bkg_eta[0]*self.bkg_eta[3])
        self.bkg_n[1] = np.sqrt(self.bkg_eta[1]*self.bkg_eta[4])
        self.bkg_n[2] = np.sqrt(self.bkg_eta[2]*self.bkg_eta[5])

        self.n_max[0] = np.sqrt(eta[0]*eta[3])
        self.n_max[1] = np.sqrt(eta[1]*eta[4])
        self.n_max[2] = np.sqrt(eta[2]*eta[5])

        return eta

    def _gaussian_rip(self,x,y,z,t=0):
        grid = [x,y,z]
        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')
        _r2 = np.zeros( [ 6,x.shape[0],y.shape[1],z.shape[2]], order='F')
        _rt = np.zeros( [ 6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        dims = dim[self.dim]

        if len(self.dim)==1:
            p = dims[0]
            for i in range(0,6):
                _temp1 = grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t
                _r2[i] = (_temp1/self.delta_sigma[p,i])**2
                _rt[i] = (_temp1*self.delta_velocity[p,i])/(self.delta_sigma[p,i]**2)

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t
                _temp2 = grid[q] - self.offset[q,i] - self.delta_velocity[q,i]*t
                _r2[i] = (_temp1/self.delta_sigma[p,i])**2 + (_temp2/self.delta_sigma[q,i])**2
                _rt[i] = (_temp1*self.delta_velocity[p,i])/(self.delta_sigma[p,i]**2) + \
                    (_temp2*self.delta_velocity[q,i])/(self.delta_sigma[q,i]**2)

        if len(self.dim)==3:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t
                _temp2 = grid[q] - self.offset[q,i] - self.delta_velocity[q,i]*t
                _temp3 = grid[r] - self.offset[r,i] - self.delta_velocity[r,i]*t
                _r2[i] = (_temp1/self.delta_sigma[p,i])**2 + (_temp2/self.delta_sigma[q,i])**2 + (_temp3/self.delta_sigma[r,i])**2
                _rt[i] = (_temp1*self.delta_velocity[p,i])/(self.delta_sigma[p,i]**2) + \
                    (_temp2*self.delta_velocity[q,i])/(self.delta_sigma[q,i]**2) + \
                    (_temp3*self.delta_velocity[r,i])/(self.delta_sigma[r,i]**2)
        _r2 = np.exp(-_r2)
        _rt = 2.0*_r2*_rt

        for i in range(0,6):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+6] = self.delta_eta[i]*_rt[i]

        return eta

    def _gaussianf(self,x,y,z):

        u = x - (self.delta_corner[0] + self.delta_smooth_length/2.0)
        v = y - (self.delta_corner[1] + self.delta_smooth_width/2.0)
        v = z - (self.delta_corner[2] + self.delta_smooth_height/2.0)

        r2 = u**2/self.delta_length**2 + v**2/self.delta_width**2 + w**2/self.delta_length**2
        
        g = np.exp(-r2)
        
        return g

    def _gaussian(self,x,y,z,t=0):
        grid = [x,y,z]
        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')
        _r2 = np.zeros( [ 6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        dims = dim[self.dim]

        if len(self.dim)==1:
            p = dims[0]
            for i in range(0,6): 
                _r2[i] = ((grid[p] - self.offset[p,i])/self.delta_sigma[p,i])**2

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i])/self.delta_sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i])/self.delta_sigma[q,i])**2
                _r2[i] = _temp1 + _temp2

        if len(self.dim)==3:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i])/self.delta_sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i])/self.delta_sigma[q,i])**2
                _temp3 = ((grid[r] - self.offset[r,i])/self.delta_sigma[r,i])**2
                _r2[i] = _temp1 + _temp2 + _temp3

        _r2 = np.exp(-_r2)

        for i in range(0,6): eta[i] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]

        return eta

    def _tanh_rip(self,x,y,z,t=0):
        grid = [x,y,z]
        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')
        _r2 = np.zeros( [ 6,x.shape[0],y.shape[1],z.shape[2]], order='F')
        _rt = _r2

        dims = dim[self.dim]

        if len(self.dim)==1:
            p = dims[0]
            for i in range(0,6):
                _r2[i] = ((grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t)/self.delta_sigma[p,i])
                _rt[i] = self.delta_velocity[p,i]/(2.0*self.delta_sigma[p,i]**2)

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t)/self.delta_sigma[p,i])
                _temp2 = ((grid[q] - self.offset[q,i] - self.delta_velocity[q,i]*t)/self.delta_sigma[q,i])
                _rt[i] = self.delta_velocity[p,i]/(2.0*self.delta_sigma[p,i]) + \
                    self.delta_velocity[q,i]/(2.0*self.delta_sigma[q,i])
                _r2[i] = _temp1 + _temp2

        if len(self.dim)==2:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.delta_velocity[p,i]*t)/self.delta_sigma[p,i])
                _temp2 = ((grid[q] - self.offset[q,i] - self.delta_velocity[q,i]*t)/self.delta_sigma[q,i])
                _temp3 = ((grid[r] - self.offset[r,i] - self.delta_velocity[r,i]*t)/self.delta_sigma[r,i])
                _rt[i] = self.delta_velocity[p,i]/(2.0*self.delta_sigma[p,i]**2) + \
                    self.delta_velocity[q,i]/(2.0*self.delta_sigma[q,i]) + \
                    self.delta_velocity[r,i]/(2.0*self.delta_sigma[r,i])
                _r2[i] = _temp1 + _temp2 + _temp3

        _r2 = 1.0 + np.tanh(_r2)
        _rt = _rt/(np.cosh(r_2)**2)

        for i in range(0,6):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+6] = self.delta_eta[i]*_rt[i]

        return eta

    def _homogeneous(self,x,y,z,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        for i in range(0,6): eta[i] = self.bkg_eta[i]

        return eta

    def _single_fiber(self,x,y,z,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        xi,yi,zi = self.fiber_corner

        spany = (y>=yi)*(y<=(yi + self.fiber_width ))
        spanx = (x>=xi)*(y<=(xi + self.fiber_length))
        spanz = (z>=zi)*(z<=(zi + self.fiber_height))
        
        span  = spanx*spany*spanz

        for i in range(0,6): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span

        return eta

    def _double_fiber(self,x,y,z,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        y1i,y2i = self.fiber_corner[:,1]
        xi1,x2i = self.fiber_corner[:,0]
        zi1,z2i = self.fiber_corner[:,2]

        span1y = ((y>=y1i)*(y<=(y1i+self.fiber_width[0])))
        span2y = ((y>=y2i)*(y<=(y2i+self.fiber_width[1])))

        span1x = ((x>=x1i)*(x<=(x1i+self.fiber_length[0])))
        span2x = ((x>=x2i)*(x<=(x2i+self.fiber_length[1])))

        span1z = ((z>=z1i)*(z<=(z1i+self.fiber_height[0])))
        span2z = ((z>=z2i)*(z<=(z2i+self.fiber_height[1])))

        span1 = span1y*span1x*span1z
        span2 = span2y*span2x*span2z

        for i in range(0,6): eta[i] = self.bkg_eta[i] + self.fiber_eta[0,i]*span1 + self.fiber_eta[1,i]*span2

        return eta

    def _oscillate(self,x,y,z,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        xi,yi,zi = self.fiber_corner

        span = ((y>=yi)*(y<=(yi+self.fiber_width)))*((x>=xi)*(x<=(xi+self.fiber_length)))*((z>=zi)*(z<=(zi+self.fiber_height)))

        xid,yid,zid = self.delta_corner
        
        spand = ((y>=yid)*(y<=(yid+self.delta_width)))*((x>=xid)*(x<=(xid+self.delta_length)))*((z>=zid)*(z<=(zid+self.delta_height)))

        w,dw = self._get_vibrate(span=spand,t=t)

        for i in range(0,6): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span + self.delta_eta[i]*w

        return eta

    def __init__(self,normalized=True,shape='homogeneous',metal=False):
        self.num_aux = 12
        self.num_dim = 3
        self.options = {}
        super(Material3D,self).__init__(normalized,shape)
        self.dim = 'x'
        self.metal = metal
        self._dx = 1.0
        self._dy = 1.0
        self._dz = 1.0
