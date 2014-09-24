import numpy as np

def user_source():
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
        if self.shape=='off':
            grid = state.grid
            if state.num_dim==1:
                x = grid.x.centers
                waveshape = self.shape_function(-(x-self.offset)**2/(self.pulse_width**2))
                state.q[0,:] = self._material.zo*waveshape
                state.q[1,:] = waveshape
            elif state.num_dim==2:
                grid = state.grid
                x,y = grid.c_centers
                waveshape = self.shape_function(-(x-self.offset[1])**2/(self.pulse_width**2))
                state.q[0,:,:] = 0.0
                state.q[1,:,:] = self._material.zo*waveshape
                state.q[2,:,:] = waveshape
        else:
            if state.num_dim==1:
                state.q[:,:] = 0.0
            elif state.num_dim==2:
                state.q[:,:,:] = 0.0
            elif state.num_dim==3:
                state.q[:,:,:,:] = 0.0

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
            self.dim = 'x'

        if self.shape=='harmonic pulse':
            self.pulse_width = self.wavelength
            self.harmonic_function = np.sin
            self.shape_function = np.exp
            self.function = self._harmonic_pulse
            self.dim = 'x'

        if self.shape=='off':
            self.pulse_width = self.wavelength
            self.shape_function = np.exp

        if self.transversal_shape=='plane':
            self.transversal_function = lambda y: 1.0
        
        if self.transversal_shape=='gauss':
            self.transversal_function = lambda y: np.exp(-(y - self.transversal_offset)**2/self.transversal_width**2)
        
        if self.transversal_shape=='cosine':
            self.transversal_function = lambda y: np.cos((y-self.transversal_offset)*np.pi/(self.transversal_width))

        return

    def _plane(self,x,y,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        harmonic = self.transversal_function(y)*self.harmonic_function(self.k[0]*(x-self.offset[1]) - self.omega*t)

        wave[0,:,:] = self.amplitude[0]*harmonic
        wave[1,:,:] = self.amplitude[1]*harmonic
        wave[2,:,:] = self.amplitude[2]*harmonic

        return wave

    def _pulse(self,x,y,t=0):

        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')
        
        shapex = self.shape_function(-(x - (self.offset[1] + self.v*t))**2/self.pulse_width**2)

        shapey = self.transversal_function(y)

        shape = shapey*shapex
        
        wave[0,:,:] = self.amplitude[0]*shape
        wave[1,:,:] = self.amplitude[1]*shape
        wave[2,:,:] = self.amplitude[2]*shape

        return wave

    def _harmonic_pulse(self,x,y,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')
        harmonic = self.harmonic_function(self.k[0]*(x-self.offset[1]) - self.omega*t)

        shape = self.transversal_function(y)*self.shape_function(-(x - (self.offset[1] + self.v*t))**2/self.pulse_width**2)
        shape = shape*harmonic

        wave[0,:,:] = self.amplitude[0]*shape
        wave[1,:,:] = self.amplitude[1]*shape
        wave[2,:,:] = self.amplitude[2]*shape

        return wave

    def _off(self,x,t=0):
        wave = np.zeros( [3,x.shape[0],y.shape[1]], order='F')

        return wave

    def __init__(self,material,shape='plane',**kwargs):
        self._set_f_w(material,kwargs)
        self.options = {}
        self.k = np.asarray([2.0*np.pi/self.wavelength,0.0])
        self.v = material.co/material.bkg_n.max()
        self.amplitude = np.asarray([0.0,material.zo,1.0])
        self.offset = np.zeros([3])
        self.transversal_shape = 'plane'
        self.transversal_offset = 0.0
        self.transversal_width = 0.0
        self.transversal_function = None
        self.shape = shape
        self.custom = False
        self.function = None
        self._material = material
