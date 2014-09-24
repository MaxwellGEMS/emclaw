import numpy as np

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

class Material:

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
                print "%s = %s" % (attr, getattr(self, attr))

    def _dump(self,obj):
        for attr in sorted(dir(obj)):
            try:
                print "%s = %s" % (attr, getattr(obj, attr))
            except:
                pass

    def set_moving_gauss(self):
        self.shape = 'moving gauss'
        self.setup()
        return

    def set_moving_tanh(self):
        self.shape = 'moving tanh'
        self.setup()
        return

    def set_gaussian(self):
        self.shape = 'gaussian'
        self.setup()
        return

    def set_custom(self):
        self.custom = True
        self.shape = 'custom'
        self.setup()
        return

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
        elif state.num_dim==2:
            x,y = grid._c_centers
            state.aux = self.function(x,y,t)
        elif state.num_dim==3:
            x,y,z = grid._c_centers
            state.aux = self.function(x,y,t)
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
        grid.compute_c_centers()
        t = state.t
        if state.num_dim==1:
            x = grid.x.centers
            state.aux = self.function(x,t)
        elif state.num_dim==2:
            x,y = grid._c_centers
            state.aux = self.function(x,y,t)
        elif state.num_dim==3:
            x,y,z = grid._c_centers
            state.aux = self.function(x,y,t)
        return state

    def __init__(self):
        self.normalized_vacuum = True
        self.shape = None
        self.custom = False
        
class Material1D(Material):
    def setup(self,options={}):
        self._unpack_options(options=options)
        self._set_vacuum()
        self.bkg_n = np.sqrt(self.bkg_er*self.bkg_mr)
        self.n_max = self.bkg_n
        temp_flag = False
        if self.custom:
            self.custom_function = user_material

        if self.shape=='homogeneous':
            self.function = self._homogeneous
        
        if self.shape.startswith('moving'):
            self.velocity_e = 0.59
            self.velocity_m = self.velocity_e
            self.offset_e   = 10.0
            self.offset_m   = self.offset_e
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
            self.sigma_m            = self.sigma_e
            self.relative_amplitude = 0.1
            self.delta_n     = self.relative_amplitude*(self.bkg_n)
            self.delta_e     = self.delta_n
            self.delta_m     = self.delta_n
            self.em_equal = True

            if not self._moving:
                self.function = self._gaussian
            
            self._rip_precalc = False

        if self.shape.find('fiber vibrate'):
            self.delta_length = 5.0
            self.delta_corner = 5.0
            self.delta_eta    = np.ones([2])
            self.delta_smooth = False
            self.delta_omega  = 2.0*np.pi
            self.delta_function = np.sin
            self.delta_smooth_function = self._gaussianf
            self.function = self._oscillate_fiber
            self.update_at_each_stage = True

        if self.nonlinear:
            self.chi2_e = 0.0
            self.chi3_e = 0.0
            self.chi2_m = 0.0
            self.chi3_m = 0.0

        return

    def _gaussian_rip(self,x,t):
                
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e = x - self.velocity_e*t - self.offset_e
        u_x_m = x - self.velocity_m*t - self.offset_m
        u_e = (u_x_e/self.sigma_e)**2
        u_m = (u_x_m/self.sigma_m)**2
        u_e_t = 2.0*((self.velocity_e*u_x_e)/(self.sigma_e**2))
        u_m_t = 2.0*((self.velocity_m*u_x_m)/(self.sigma_m**2))

        eta[0,:] = self.delta_e*np.exp(-u_e) + self.bkg_er
        eta[1,:] = self.delta_m*np.exp(-u_m) + self.bkg_mr
        eta[2,:] = u_e_t*self.delta_e*np.exp(-u_e)
        eta[3,:] = u_m_t*self.delta_m*np.exp(-u_m)

        return eta

    def _gaussian(self,x,t=0):
        
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e = x - self.offset_e
        u_x_m = x - self.offset_m
        u_e = (u_x_e/selfsigma_e)**2
        u_m = (u_x_m/selfsigma_m)**2
        eta[0,:] = self.delta_e*np.exp(-u_e) + self.bkg_er
        eta[1,:] = self.delta_m*np.exp(-u_m) + self.bkg_mr    

        return eta

    def _gaussianf(self,x,y):

        u = x - self.delta_corner + self.delta_length/2.0

        r2 = u**2/self.delta_length**2
        
        g = np.exp(-r2)
        
        return g

    def _oscillate_fiber(self,x,t=0):

        eta = np.zeros( [4,x.shape[0],y.shape[1]], order='F')

        xid = self.delta_corner
        
        spand = ((x>=xid)*(x<=(xid+self.delta_length)))

        if not self.delta_smooth:
            w = self.delta_function(self.delta_omega*t)*spand
        else:
            w = self.delta_function(self.delta_omega*t)*self.delta_smooth_function(x)*spand

        for i in range(0,2): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span + self.delta_eta[i]*w

        return eta

    def _tanh_rip(self,x,t):
        
        eta = np.zeros( [4,len(x)], order='F')

        u_x_e = (x - self.velocity_e*t - self.offset_e)/self.sigma_e
        u_x_m = (x - self.velocity_m*t - self.offset_m)/self.sigma_m

        eta[0,:] = (self.delta_e/2.0)*(1.0 + np.tanh(u_x_e)) + self.bkg_er
        eta[1,:] = (self.delta_m/2.0)*(1.0 + np.tanh(u_x_m)) + self.bkg_mr
        eta[2,:] = -(self.delta_e*self.velocity_e/(2.0*self.sigma_e))/(np.cosh(u_x_e)**2)
        eta[3,:] = -(self.delta_m*self.velocity_m/(2.0*self.sigma_m))/(np.cosh(u_x_m)**2)

        return eta

    def _homogeneous(self,x,t=0):

        eta = np.zeros( [4,len(x)], order='F')

        eta[0,:] = self.bkg_er
        eta[1,:] = self.bkg_mr

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
        self.normalized_vacuum = normalized
        self.bkg_er = 1.0
        self.bkg_mr = 1.0
        self.shape = shape
        self.options = {}
        self.nonlinear = True
        self.custom = False
        self._moving = False

class Material2D(Material):
    def setup(self,options={}):
        self._unpack_options(options=options)
        self._set_vacuum()
        self.bkg_n = np.ones([2])
        self.n_max = np.ones([2])

        temp_flag = False
        if self.custom:
            self.custom_function = user_material

        if self.shape=='homogeneous':
            self.function = self._homogeneous
        
        if self.shape.startswith('moving'):
            self.velocity = np.append(0.59*np.ones([1,3]),np.zeros([1,3]),axis=0)
            self.offset   = np.append(10.0*np.ones([1,3]),np.zeros([1,3]),axis=0)
            self._moving    = True
            self.update_at_each_stage = True

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
            self.delta_sigma = 5.0*np.ones([2,3])
            self.relative_amplitude = 0.1*np.ones([3])
            self.delta_eta = self.relative_amplitude*self.bkg_eta
            self.em_equal = True

            if not self._moving:
                self.function = self._gaussian
            
            self._rip_precalc = False

        if self.shape.startswith('fiber'):
            self.fiber_eta = np.ones([3])

        if self.shape=='fiber single':
            self.fiber_corner = [-5.0,0.0]
            self.fiber_width  = 5.0
            self.fiber_length = 100.0
            self.function = self._single_fiber

        if self.shape=='fiber double':
            self.fiber_eta    = np.ones([2,3])
            self.fiber_corner = np.zeros([2,2])
            self.fiber_width  = np.ones([2])
            self.fiber_length = 100.0*np.ones([2])
            self.function = self._double_fiber

        if self.shape=='fiber vibrate':
            self.fiber_corner = [-5.0,0.0]
            self.fiber_width  = 5.0
            self.fiber_length = 100.0
            self.delta_width  = 5.0
            self.delta_length = 5.0
            self.delta_corner = [5.0,0.0]
            self.delta_eta    = np.ones([3])
            self.delta_smooth = False
            self.delta_omega  = 2.0*np.pi
            self.delta_function = np.cos
            self.delta_smooth_function = self._gaussianf
            self.function = self._oscillate_fiber
            self.update_at_each_stage = True
        
        if self.nonlinear:
            self.chi2 = np.zeros( [3], order='F')
            self.chi3 = np.zeros( [3], order='F')

        if self.metal:
            self.metal_corners = []

        return

    def set_fiber_single(self):
        self.shape = 'fiber single'
        self.setup()
        return

    def set_fiber_double(self):
        self.shape = 'fiber double'
        self.setup()
        return

    def set_fiber_single_sine(self):
        self.shape = 'fiber single sine'
        self.setup()
        return

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
        _rt = _r2

        u = v = False

        if self.dim=='x': u = True
        
        if self.dim=='y': v =  True
        
        if self.dim=='xy': u = v = True

        if u:
            for i in range(0,3): 
                _r2[i] = ((x - self.offset[0,i] - self.velocity[0,i]*t)/self.delta_sigma[0,i])**2
                _rt[i] = 2.0*(self.velocity[0,i]*_r2[i])/(self.delta_sigma[0,i]**2)

        if v:
            for i in range(0,3):
                _temp  = ((y - self.offset[1,i] - self.velocity[1,i]*t)/self.delta_sigma[1,i])**2
                _rt[i] = _rt[i] + 2.0*(self.velocity[1,i]*_temp)/(self.delta_sigma[1,i]**2)
                _r2[i] = _r2[i] + _temp
        _r2 = np.exp(-_r2)
        _rt = _r2*_rt

        for i in range(0,3):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+3] = self.delta_eta[i]*_rt[i]

        return eta

    def _gaussianf(self,x,y):

        u = x - self.delta_corner[0] + self.delta_length/2.0
        v = y - self.delta_corner[1] + self.delta_width/2.0

        r2 = u**2/self.delta_length**2 + v**2/self.delta_width**2
        
        g = np.exp(-r2)
        
        return g

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
            for i in range(0,3): _r2[i] = (x - self.velocity[0,i] - self.offset[0,i])/self.delta_sigma[0,i]
            _rt = self.velocity[0,:]/(2.0*self.delta_sigma[0,:])
        if v:
            for i in range(0,3): _r2[i] = _r2[i] + (y - self.velocity[1,i] - self.offset[1,i])/self.delta_sigma[1,i]
            _rt = _rt + self.velocity[1,:]/(2.0*self.delta_sigma[1,:])

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

    def _oscillate_fiber(self,x,y,t=0):

        eta = np.zeros( [6,x.shape[0],y.shape[1]], order='F')

        xi,yi = self.fiber_corner

        span = ((y>=yi)*(y<=(yi+self.fiber_width)))*((x>=xi)*(x<=(xi+self.fiber_length)))

        xid,yid = self.delta_corner
        
        spand = ((y>=yid)*(y<=(yid+self.delta_width)))*((x>=xid)*(x<=(xid+self.delta_length)))

        if not self.delta_smooth:
            w = self.delta_function(self.delta_omega*t)*spand
        else:
            w = self.delta_function(self.delta_omega*t)*self.delta_smooth_function(x,y)*spand

        for i in range(0,3): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span + self.delta_eta[i]*w

        return eta

    def __init__(self,normalized=True,shape='homogeneous',metal=False):
        self.normalized_vacuum = normalized
        self.bkg_eta = np.ones([3])
        self.shape = shape
        self.options = {}
        self.nonlinear = True
        self.custom = False
        self._moving = False
        self.dim = 'x'
        self.v = 1.0
        self.update_at_each_stage = False
        self.metal = metal

class Material3D(Material):
    def setup(self,options={}):
        self._unpack_options(options=options)
        self._set_vacuum()
        self.bkg_n = np.ones([3])
        self.n_max = np.ones([3])

        temp_flag = False
        if self.custom:
            self.custom_function = user_material

        if self.shape=='homogeneous':
            self.function = self._homogeneous
        
        if self.shape.startswith('moving'):
            self.velocity = np.append(0.59*np.ones([1,3]),np.zeros([1,3]),np.zeros([1,3]),axis=0)
            self.offset   = np.append(10.0*np.ones([1,3]),np.zeros([1,3]),np.zeros([1,3]),axis=0)
            self._moving  = True
            self.update_at_each_stage = True

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
            self.sigma = 5.0*np.ones([3,3])
            self.relative_amplitude = 0.1*np.ones([6])
            self.delta_eta = self.relative_amplitude*self.eta
            self.em_equal = True

            if not self._moving:
                self.function = self._gaussian
            
            self._rip_precalc = False

        if self.shape.startswith('fiber'):
            self.fiber_eta = np.ones([6])

        if self.shape=='fiber single':
            self.fiber_corner = [-5.0,0.0,0.0]
            self.fiber_width  = 5.0
            self.fiber_height = 5.0
            self.fiber_length = 100.0
            self.function = self._single_fiber

        if self.shape=='fiber double':
            self.fiber_corner = np.append([5.0,0,0],[0,0,0],axis=1).reshape(2,3)
            self.fiber_width  = np.ones([2])
            self.fiber_height = np.ones([2])
            self.fiber_length = 100.0*np.ones([2])
            self.function = self._double_fiber

        if self.shape.find('fiber vibrate'):
            self.fiber_corner = [-5.0,0.0,0.0]
            self.fiber_width  = 5.0
            self.fiber_height = 5.0
            self.fiber_length = 100.0
            self.delta_width  = 5.0
            self.delta_height = 5.0
            self.delta_length = 100.0
            self.delta_corner = [5.0,0.0,0.0]
            self.delta_eta    = np.ones([6])
            self.delta_smooth = False
            self.delta_omega  = 2.0*np.pi
            self.delta_function = np.cos
            self.delta_smooth_function = self._gaussianf
            self.function = self._oscillate_fiber
            self.update_at_each_stage = True
        
        if self.nonlinear:
            self.chi2 = np.zeros( [6], order='F')
            self.chi3 = np.zeros( [6], order='F')

        if self.metal:
            self.metal_corners = []

        return

    def set_fiber_single(self):
        self.shape = 'fiber single'
        self.setup()
        return

    def set_fiber_double(self):
        self.shape = 'fiber double'
        self.setup()
        return

    def set_fiber_single_sine(self):
        self.shape = 'fiber single sine'
        self.setup()
        return

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
        _rt = _r2

        dims = dim[self.dim]

        if len(self.dim)==1:
            p = dims[0]
            for i in range(0,6): 
                _r2[i] = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2
                _rt[i] = 2.0*(self.velocity[p,i]*_r2[i])/(self.sigma[p,i]**2)

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])**2
                _rt[i] = 2.0*(self.velocity[p,i]*_temp1)/(self.sigma[p,i]**2) + 2.0*(self.velocity[q,i]*_temp2)/(self.sigma[q,i]**2)
                _r2[i] = _temp1 + _temp2

        if len(self.dim)==3:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])**2
                _temp3 = ((grid[r] - self.offset[r,i] - self.velocity[r,i]*t)/self.sigma[r,i])**2
                _rt[i] = 2.0*(self.velocity[p,i]*_temp1)/(self.sigma[p,i]**2) + 2.0*(self.velocity[q,i]*_temp2)/(self.sigma[q,i]**2) + 2.0*(self.velocity[r,i]*_temp3)/(self.sigma[r,i]**2)
                _r2[i] = _temp1 + _temp2 + _temp3


        _r2 = np.exp(-_r2)
        _rt = _r2*_rt

        for i in range(0,6):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+6] = self.delta_eta[i]*_rt[i]

        return eta

    def _gaussianf(self,x,y,z):

        u = x - self.delta_corner[0] + self.delta_width[0]/2.0
        v = y - self.delta_corner[1] + self.delta_width[1]/2.0
        v = z - self.delta_corner[2] + self.delta_width[2]/2.0

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
                _r2[i] = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])**2
                _r2[i] = _temp1 + _temp2

        if len(self.dim)==3:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])**2
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])**2
                _temp3 = ((grid[r] - self.offset[r,i] - self.velocity[r,i]*t)/self.sigma[r,i])**2
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
                _r2[i] = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])
                _rt[i] = self.velocity[p,i]/(2.0*self.sigma[p,i]**2)

        if len(self.dim)==2:
            p,q = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])
                _rt[i] = self.velocity[p,i]/(2.0*self.sigma[p,i]) + self.velocity[q,i]/(2.0*self.sigma[q,i])
                _r2[i] = _temp1 + _temp2

        if len(self.dim)==2:
            p,q,r = dims
            for i in range(0,6):
                _temp1 = ((grid[p] - self.offset[p,i] - self.velocity[p,i]*t)/self.sigma[p,i])
                _temp2 = ((grid[q] - self.offset[q,i] - self.velocity[q,i]*t)/self.sigma[q,i])
                _temp3 = ((grid[r] - self.offset[r,i] - self.velocity[r,i]*t)/self.sigma[r,i])
                _rt[i] = self.velocity[p,i]/(2.0*self.sigma[p,i]**2) + self.velocity[q,i]/(2.0*self.sigma[q,i]) + self.velocity[r,i]/(2.0*self.sigma[r,i])
                _r2[i] = _temp1 + _temp2 + _temp3

        _r2 = 1.0 + np.tanh(_r2)
        _rt = _rt/(np.cosh(r_2)**2)

        for i in range(0,6):
            eta[i  ] = self.delta_eta[i]*_r2[i] + self.bkg_eta[i]
            eta[i+6] = self.delta_eta[i]*_rt[i]

        return eta

    def _homogeneous(self,x,y,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        for i in range(0,6): eta[i] = self.bkg_eta[i]

        return eta

    def _single_fiber(self,x,y,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        xi,yi,zi = self.fiber_corner

        spany = (y>=yi)*(y<=(yi + self.fiber_width ))
        spanx = (x>=xi)*(y<=(xi + self.fiber_length))
        spanz = (z>=zi)*(z<=(zi + self.fiber_height))
        
        span  = spanx*spany*spanz

        for i in range(0,6): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span

        return eta

    def _double_fiber(self,x,y,t=0):

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

    def _oscillate_fiber(self,x,y,z,t=0):

        eta = np.zeros( [12,x.shape[0],y.shape[1],z.shape[2]], order='F')

        xi,yi,zi = self.fiber_corner

        span = ((y>=yi)*(y<=(yi+self.fiber_width)))*((x>=xi)*(x<=(xi+self.fiber_length)))*((z>=zi)*(z<=(zi+self.fiber_height)))

        eta[0,:,:] = self.bkg_eta[0] + self.fiber_eta[0]*span
        eta[1,:,:] = self.bkg_eta[1] + self.fiber_eta[1]*span
        eta[2,:,:] = self.bkg_eta[2] + self.fiber_eta[2]*span

        xid,yid,zid = self.delta_corner
        
        spand = ((y>=yid)*(y<=(yid+self.delta_width)))*((x>=xid)*(x<=(xid+self.delta_length)))*((z>=zid)*(z<=(zid+self.delta_height)))

        if not self.delta_smooth:
            w = self.delta_function(self.delta_omega*t)*spand
        else:
            w = self.delta_function(self.delta_omega*t)*self.delta_smooth_function(x,y,z)*spand

        for i in range(0,6): eta[i] = self.bkg_eta[i] + self.fiber_eta[i]*span + self.delta_eta[i]*w

        return eta

    def __init__(self,normalized=True,shape='homogeneous',metal=False):
        self.normalized_vacuum = normalized
        self.bkg_eta = np.ones([6])
        self.shape = shape
        self.options = {}
        self.nonlinear = True
        self.custom = False
        self._moving = False
        self.dim = 'x'
        self.v = 1.0
        self.update_at_each_stage = False
        self.metal = metal
        