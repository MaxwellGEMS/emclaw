import numpy as np


heading = {}
heading['xyz'] = [0,1,2]
heading['xy' ] = [0,1]
heading['xz' ] = [0,2]
heading['yz' ] = [1,2]
heading['x'  ] = [0]
heading['y'  ] = [1]
heading['z'  ] = [2]

polarization = {}
polarization['xyz'] = [0,1,2]
polarization['xy' ] = [1,0]
polarization['xz' ] = [2,0]
polarization['yz' ] = [2,1]
polarization['x'  ] = [1]
polarization['y'  ] = [0]
polarization['z'  ] = [0]

def user_source():
    pass

def user_transversal_source():
    pass

class Sources:
    def _unpack_options(self,options={}):
        # first add the options passed to the function
        for key in options:
            setattr(self,key,options[key])

        # unpack options from self.options{}
        for key in self.options:
            setattr(self,key,self.options[key])

    def _set_f_w(self,material,dictin):
        
        setattr(self,'c',material.co)

        for key,value in dictin.iteritems():
            setattr(self,key,value)

        if hasattr(self,'wavelength'):
            setattr(self,'omega',2.0*np.pi*material.co/self.wavelength)
        elif hasattr(self,'frequency'):
            setattr(self,'wavelength',2.0*np.pi*material.co/self.frequency)
        else:
            msg = 'You must define either wavelength or frequency'
            print msg

    def set_plane(self):
        self.shape = 'plane'
        self.setup()
        return

    def set_pulse(self):
        self.shape = 'pulse'
        self.setup()
        return

    def set_harmonic_pulse(self):
        self.shape = 'harmonic pulse'
        self.setup()
        return

    def set_custom(self):
        self.custom = True
        self.setup()

    def set_off(self):
        self.shape = 'off'
        self.setup()
        return

    def init(self,state):
        state.q.fill(0.0)

        if self.shape=='off':

            dimh = heading[self.heading]
            dimp = polarization[self.heading]

            grid = state.grid

            if state.num_dim==1:
                x = grid.x.centers
                waveshape = np.exp(-(x-self.offset)**2/(self.pulse_width**2))
                state.q[0,:] = self._material.zo*waveshape
                state.q[1,:] = waveshape
            if state.num_dim>=2:
                _r = 0.0
                for i in range(len(self.heading)):
                    _r = _r + (grid.c_centers[dimh[i]] - self.offset[dimp[i]])**2/self.pulse_width**2

                waveshape = np.exp(-_r)

                if len(self.heading)==1:
                    state.q[dimp[0],:,:] = ((-1.0)**dimh[0])*self._material.zo*waveshape
                
                if state.num_dim==2: p=2
                if state.num_dim==3: p=5
                
                state.q[p,:,:] = waveshape

        if self.shape=='farago':
            x = state.grid.x.centers
            dx = state.grid.x.delta
            ddx = dx/2.0
            state.q[0,:] = (2.0*np.sin(ddx)*np.sin(x))/dx      
            #state.q[1,:] = #-1.0*(np.cos(x)*np.sin(ddx))/dx

        return state

    def scattering_bc(self,state,dim,t,qbc,auxbc,num_ghost):
        grid = state.grid
        grid.compute_c_centers_with_ghost(num_ghost,recompute=True)
        t = state.t
        if state.num_dim==1:
            x = grid.x.centers_with_ghost[:num_ghost]
            qbc[:,:num_ghost] = self.function(x,t)

        if state.num_dim==2:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            if dim.name==state.grid.dimensions[0].name:
                x = x[:num_ghost,:]
                y = y[:num_ghost,:]
                qbc[:,:num_ghost,:] = self.function(x,y,t)
            else:
                x = x[:,:num_ghost]
                y = y[:,:num_ghost]
                qbc[:,:,:num_ghost] = self.function(x,y,t)

        if state.num_dim==3:
            x = grid._c_centers_with_ghost[0]
            y = grid._c_centers_with_ghost[1]
            z = grid._c_centers_with_ghost[2]
            if dim.name==state.grid.dimensions[0].name:
                x = x[:num_ghost,:,:]
                y = y[:num_ghost,:,:]
                z = z[:num_ghost,:,:]
                qbc[:,:num_ghost,:,:] = self.function(x,y,z,t)
            if dim.name==state.grid.dimensions[1].name:
                x = x[:,:num_ghost,:]
                y = y[:,:num_ghost,:]
                z = z[:,:num_ghost,:]
                qbc[:,:,:num_ghost,:] = self.function(x,y,z,t)
            if dim.name==state.grid.dimensions[2].name:
                x = x[:,:,:num_ghost]
                y = y[:,:,:num_ghost]
                z = z[:,:,:num_ghost]
                qbc[:,:,:,:num_ghost] = self.function(x,y,z,t)

        return qbc

    def dump(self):
        for attr in sorted(dir(self)):
            if not attr.startswith('_'):
                print "%s = %s" % (attr, getattr(self, attr))

    def __init__(self):
        self.shape      = None
        self.custom     = False
        self.custom_func = user_source

class Source1D(Sources):

    def setup(self,options={}):
        self._unpack_options(options=options)

        if self.shape=='plane':
            self.harmonic_function = np.sin
            self.function = self._plane

        if self.shape=='pulse':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp
            self.function = self._pulse
            self.averaged = True
            self._dx = 1.0
            self._cp = np.sqrt(np.pi)

        if self.shape=='harmonic pulse':
            self.pulse_width = self.wavelength
            self.harmonic_function = np.sin
            self.shape_function = np.exp
            self.function = self._harmonic_pulse

        if self.shape=='off':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp

        return

    def _plane(self,x,t):
        wave = np.zeros( [2,len(x)], order='F')

        harmonic = self.harmonic_function(self.k*(x-self.offset) - self.frequency*t)

        wave[0,:] = self.Ey*harmonic
        wave[1,:] = self.Hz*harmonic

        return wave

    def _pulse(self,x,t):
        wave = np.zeros( [2,len(x)], order='F')
        if self.averaged:
            from scipy.special import erf
            ddx = self._dx/2.0
            arg0 = self.offset + self.v*t - x
            arg1 = (ddx + arg0)/self.pulse_width
            arg2 = (ddx - arg0)/self.pulse_width
            pulseshape = np.sqrt(np.pi)*self.pulse_width*(erf(arg1)+erf(arg2))/(2.0*self._dx)
        else:
            pulseshape = self.shape_function(-(x - (self.offset + self.v*t))**2/self.pulse_width**2)
        wave[0,:] = self.Ey*pulseshape
        wave[1,:] = self.Hz*pulseshape

        return wave

    def _harmonic_pulse(self,x,t):

        wave = np.zeros( [2,len(x)], order='F')
        harmonic = self.harmonic_function(self.k*(x-self.offset) - self.frequency*t)
        pulseshape = self.shape_function(-(x - (self.offset + self.v*t))**2/self.pulse_width**2)
        wave[0,:] = self.Ey*harmonic*pulseshape
        wave[1,:] = self.Hz*harmonic*pulseshape

        return wave

    def _off(self,x,t=0):
        wave = np.zeros( [2,len(x)], order='F')

        return wave

    def __init__(self,material,shape='plane',**kwargs):
        self._set_f_w(material,kwargs)
        self.options = {}
        self.k = 2.0*np.pi/self.wavelength
        self.v = material.co/material.bkg_n
        self.Ey = material.zo
        self.Hz = 1.0
        self.offset = 0.0
        self.shape = shape
        self.custom = False
        self.function = None
        self.custom_func = user_source
        self._material = material
        self.dy = 1.0
        self.dx = 1.0

class Source2D(Sources):
    def setup(self,options={}):
        self._unpack_options(options=options)

        if self.shape=='custom':
            self.custom = True

        if self.custom:
            self.shape = 'custom'
            self.custom_function = user_source

        if self.shape=='plane':
            self.harmonic_function = np.sin
            self.function = self._plane

        if self.shape=='pulse':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp
            self.function = self._pulse
            self.heading = 'x'
            self.t_off = (4.0*self.pulse_width)/self.v[0]
            self.averaged = False

        if self.shape=='harmonic pulse':
            self.pulse_width = self.wavelength
            self.harmonic_function = np.sin
            self.shape_function = np.exp
            self.function = self._harmonic_pulse
            self.heading = 'x'

        if self.shape=='bessel pulse':
            self.pulse_width = self.wavelength
            self.bessel_order = 0
            self.function = self._bessel_pulse
            self.kill_after_first_zero = True

        if self.shape=='off':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp

        if self.transversal_shape=='plane':
            self.transversal_function = lambda y: 1.0

        if self.transversal_shape=='gauss':
            self.transversal_function = lambda y: self._transversal_gauss(y)

        if self.transversal_shape=='cosine':
            self.transversal_function = lambda y: self._transversal_cosine(y)

        if self.transversal_shape=='bessel':
            self.transversal_bessel_order = 0
            self.transversal_kill_after_first_zero = True
            self.transversal_function = self._transversal_bessel

        return

    def _trasversal_plane(self,u):
        shape = 1.0
        return shape

    def _transversal_cosine(self,u):
        r = (u-self.transversal_offset)/self.transversal_width
        shape = np.cos(r*np.pi)*(np.abs(r)<=0.5)
        return shape

    def _transversal_gauss(self,u):
        r = (u - self.transversal_offset)**2/self.transversal_width**2
        shape = np.exp(-r)
        return shape

    def _plane(self,x,y,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        harmonic = self.transversal_function(y)*self.harmonic_function(self.k[0]*(x-self.offset[0] - (self.omega/self.k[0])*t))

        wave[0,:,:] = self.amplitude[0]*harmonic
        wave[1,:,:] = self.amplitude[1]*harmonic
        wave[2,:,:] = self.amplitude[2]*harmonic

        return wave

    def _pulse(self,x,y,t=0):

        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        if t<=self.t_off:
            if self.averaged:
                shape = self.transversal_function(y)*self._pulse_averaged(x,y,t=t)
            else:
                shape = self.transversal_function(y)*np.exp(-(x - (self.offset[0] + self.v[0]*t))**2/self.pulse_width**2)
        else:
            shape = 0.0
        p = None
        if self.heading=='x': p=1
        if self.heading=='y': p=0
        
        if p is not None:
            wave[p,:,:] = self.amplitude[p]*shape
        
        wave[2,:,:] = self.amplitude[2]*shape

        return wave

    def _harmonic_pulse(self,x,y,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        if t<=self.t_off:
            harmonic = self.harmonic_function(self.k[0]*(x-self.offset[1] - (self.omega/self.k[0])*t))
            shape = self.transversal_function(y)*self.shape_function(-(x - (self.offset[1] + self.v[0]*t))**2/self.pulse_width**2)
            shape = shape*harmonic
        else:
            shape = 0.0

        wave[0,:,:] = self.amplitude[0]*shape
        wave[1,:,:] = self.amplitude[1]*shape
        wave[2,:,:] = self.amplitude[2]*shape

        return wave

    def _bessel_pulse(self,x,y,t=0):
        from scipy.special import jn, jn_zeros
        first_zero = jn_zeros(self.bessel_order,1)

        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        shapex = jn(self.bessel_order,(x - (self.offset[1] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))

        if self.kill_after_first_zero:
            shape_kill = np.abs((x - (self.offset[1] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))<=(first_zero[0])
            shapex = shape_kill*shapex

        shapey = self.transversal_function(y)

        shape = shapey*shapex

        wave[0,:,:] = self.amplitude[0]*shape
        wave[1,:,:] = self.amplitude[1]*shape
        wave[2,:,:] = self.amplitude[2]*shape

        return wave

    def _off(self,x,y,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        return wave

    def _transversal_bessel(self,y):
        from scipy.special import jn, jn_zeros
        first_zero = jn_zeros(self.transversal_bessel_order,1)
        shape = jn(self.transversal_bessel_order,(y-self.transversal_offset)*(first_zero[0])/(self.transversal_width/2.0))
        if self.transversal_kill_after_first_zero:
            shape_kill = np.abs((y-self.transversal_offset)*(first_zero[0])/(self.transversal_width/2.0))<=(first_zero[0])
            shape = shape_kill*shape

        return shape

    def _averaged_gauss(self,x,delta=None,s=1.0,xo=0.0,v=0.0,t=0.0):
        from scipy.special import erf
        arg = xo + v*t - x
        if delta is None:
            delta = self._dx
        ddx = delta/2.0

        erravg = (np.sqrt(np.pi)*s*(erf((ddx + arg)/s) + erf((ddx - arg)/s)))/(2.0*delta)
        
        return erravg

    def _pulse_averaged(self,x,y,t=0):
        shape = np.ones(  [x.shape[0],y.shape[1]], order='F')

        u = v = False

        if self.heading=='x': u = True
        
        if self.heading=='y': v =  True
        
        if self.heading=='xy':
            u = v = True
            self.v[1] = 0.0

        if u:
            shape = self._averaged_gauss(x,delta=self._dx,s=self.pulse_width,xo=self.offset[0],v=self.v[0],t=t)

        if v:
            shape = (shape)*self._averaged_gauss(y,delta=self._dy,s=self.pulse_width,xo=self.offset[1],v=self.v[1],t=t)

        return shape

    def __init__(self,material,shape='off',**kwargs):
        self._set_f_w(material,kwargs)
        self.options = {}
        self.k = np.asarray([2.0*np.pi/self.wavelength,0.0])
        self.v = material.co*np.asarray([1.0/material.bkg_n[0],1.0/material.bkg_n[1]])
        self.amplitude = np.asarray([0.0,material.zo,1.0])
        self.offset = np.zeros([2])
        self.transversal_shape = 'plane'
        self.transversal_offset = 0.0
        self.transversal_width = 0.0
        self.transversal_function = lambda y: 1.0
        self.shape = shape
        self.custom = False
        self.function = None
        self.heading = 'x'
        self._dx = 1.0
        self._dy = 1.0
        self._material = material

class Source3D(Sources):
    def setup(self,options={}):
        self._unpack_options(options=options)

        if self.shape=='custom':
            self.custom = True

        if self.custom:
            self.shape = 'custom'
            self.custom_function = user_source

        if self.shape=='plane':
            self.harmonic_function = np.sin
            self.function = self._plane

        if self.shape=='pulse':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp
            self.function = self._pulse
            self.heading = 'x'

        if self.shape=='harmonic pulse':
            self.pulse_width = self.wavelength
            self.harmonic_function = np.sin
            self.shape_function = np.exp
            self.function = self._harmonic_pulse
            self.heading = 'x'

        if self.shape=='bessel pulse':
            self.pulse_width = self.wavelength
            self.bessel_order = 0
            self.function = self._bessel_pulse
            self.kill_after_first_zero = True

        if self.shape=='off':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp

        if self.transversal_shape=='plane':
            self.transversal_function = lambda y,z: 1.0

        if self.transversal_shape=='gauss':
            self.transversal_function = lambda y,z: self._transversal_gauss(y,0)*self._transversal_gauss(z,1)

        if self.transversal_shape=='cosine':
            self.transversal_function = lambda y,z: self._transversal_cosine(y,0)*self._transversal_cosine(z,1)

        if self.transversal_shape=='bessel':
            self.transversal_bessel_order = 0
            self.transversal_kill_after_first_zero = True
            self.transversal_function = lambda y,z: self._transversal_bessel(y,0)*self._transversal_bessel(z,1)

        return

    def _transversal_cosine(self,u,p):
        r = (u-self.transversal_offset[p])/(self.transversal_width[p])
        shape = np.cos(r*np.pi)*(np.abs(r)<=0.5)
        return shape

    def _transversal_gauss(self,u,p):
        r = (u - self.transversal_offset[p])**2/self.transversal_width[p]**2
        shape = np.exp(-r)
        return shape

    def _transversal_bessel(self,u,p):
        from scipy.special import jn, jn_zeros
        first_zero = jn_zeros(self.transversal_bessel_order,1)
        r = (u-self.transversal_offset[p])/(self.transversal_width[p]/2.0)
        shape = jn(self.transversal_bessel_order,r*(first_zero[0]))
        if self.transversal_kill_after_first_zero:
            shape_kill = np.abs(r*(first_zero[0]))<=(first_zero[0])
            shape = shape_kill*shape
        return shape

    def _plane(self,x,y,z,t=0):
        wave = np.zeros( [6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        harmonic1 = self.transversal_function(y,z)*self.harmonic_function(self.k[0]*(x-self.offset[1]) - self.omega*t)
        harmonic2 = self.transversal_function(y,z)*self.harmonic_function(self.k[0]*(x-self.offset[5]) - self.omega*t)

        wave[1,:,:,:] = self.amplitude[1]*harmonic1
        wave[5,:,:,:] = self.amplitude[5]*harmonic2

        return wave

    def _pulse(self,x,y,z,t=0):

        wave = np.zeros( [6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        shapex1 = self.transversal_function(y,z)*self.shape_function(-(x - (self.offset[1] + self.v[0]*t))**2/self.pulse_width**2)
        shapex2 = self.transversal_function(y,z)*self.shape_function(-(x - (self.offset[5] + self.v[0]*t))**2/self.pulse_width**2)

        wave[1,:,:,:] = self.amplitude[1]*shapex1
        wave[5,:,:,:] = self.amplitude[5]*shapex2

        return wave

    def _harmonic_pulse(self,x,y,z,t=0):
        wave = np.zeros( [6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        harmonic1 = self.harmonic_function(self.k[0]*(x-self.offset[1]) - self.omega*t)
        harmonic2 = self.harmonic_function(self.k[0]*(x-self.offset[5]) - self.omega*t)

        shapex1 = self.transversal_function(y,z)*self.shape_function(-(x - (self.offset[1] + self.v[0]*t))**2/self.pulse_width**2)
        shapex2 = self.transversal_function(y,z)*self.shape_function(-(x - (self.offset[5] + self.v[0]*t))**2/self.pulse_width**2)

        wave[1,:,:,:] = self.amplitude[1]*shapex1*harmonic1
        wave[5,:,:,:] = self.amplitude[5]*shapex2*harmonic2

        return wave

    def _bessel_pulse(self,x,y,z,t=0):
        from scipy.special import jn, jn_zeros
        first_zero = jn_zeros(self.bessel_order,1)

        wave = np.zeros( [6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        shapex1 = self.transversal_function(y,z)*jn(self.bessel_order,(x - (self.offset[1] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))
        shapex2 = self.transversal_function(y,z)*jn(self.bessel_order,(x - (self.offset[5] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))

        if self.kill_after_first_zero:
            shape_kill1 = np.abs((x - (self.offset[1] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))<=(first_zero[0])
            shape_kill2 = np.abs((x - (self.offset[5] + self.v[0]*t)*(first_zero[0])/(self.pulse_width/2.0)))<=(first_zero[0])

            shapex1 = shape_kill1*shapex1
            shapex2 = shape_kill2*shapex2

        wave[1,:,:,:] = self.amplitude[1]*shapex1
        wave[5,:,:,:] = self.amplitude[5]*shapex2

        return wave

    def _off(self,x,y,z,t=0):
        wave = np.zeros( [6,x.shape[0],y.shape[1],z.shape[2]], order='F')

        return wave

    def __init__(self,material,shape='off',**kwargs):
        self._set_f_w(material,kwargs)
        self.options = {}
        self.k = np.asarray([2.0*np.pi/self.wavelength,0.0,0.0])
        self.v = material.co*np.asarray([1.0/material.bkg_n[0],1.0/material.bkg_n[1],1.0/material.bkg_n[2]])
        self.amplitude = np.asarray([0.0,material.zo,0.0,0.0,0.0,1.0])
        self.offset = np.zeros([6])
        self.transversal_shape = shape
        self.transversal_offset = [0.0,0.0]
        self.transversal_width = [0.0,0.0]
        self.transversal_function = None
        self.shape = shape
        self.custom = False
        self.function = None
        self._material = material
