_mpl_load = False
fontsize = 18
from clawpack.petclaw import Solution
from petsc4py import PETSc
from glob import glob
import os
import numpy as np
import pickle
try:
    import matplotlib
    matplotlib.use('Agg')
    matplotlib.rcParams.update({'font.size': fontsize})
    matplotlib.rcParams.update({'font.weight': 'normal'})
    matplotlib.rcParams['axes.formatter.limits'] = [0,3]
    matplotlib.rcParams['mathtext.default'] = 'sf'
    matplotlib.rcParams['axes.formatter.use_mathtext'] = True
    matplotlib.rcParams['xtick.labelsize'] = fontsize
    matplotlib.rcParams['ytick.labelsize'] = fontsize
    matplotlib.rcParams['axes.labelsize'] = fontsize
    matplotlib.rcParams['lines.linewidth'] = 1.5
    matplotlib.rcParams['lines.markersize'] = 5.0
    matplotlib.rcParams['lines.color'] = 'r'
    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import colorConverter
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    from matplotlib.lines import Line2D
    from matplotlib import pylab as plt
    from clawpack.pyclaw import Solution
    from scipy.io import loadmat,savemat
    from matplotlib.streamplot import  streamplot
    _mpl_load = True
except:
    pass

class toolbox(object):
    """docstring for VisTools"""
    def __init__(self):
        self.frame_range = None
        self.aux_range = None
        self.q_range   = None
        self.s_range   = None
        self.p_range   = None
        self.num_dim   = None
        self.num_mp    = None
        self.num_aux   = None
        self.num_eqn   = None
        self.num_cells = None
        self.num_frames = None
        self.file_format = 'petsc'
        self.outdir    = './_output'
        self._pkl_base = 'claw.pkl'
        self._ptc_base = 'claw.ptc'
        self._aux_base = 'claw_aux.ptc'
        self.xmf_base  = 'claw.meta'
        self._p_base   = 'claw_p.ptc'
        self._p_pkl    = 'claw_p.pkl'
        self._basic_info = ['names','num_cells','delta','lower','num_dim','num_eqn','num_aux']
        self._info_loaded = False
        self._grid     = False
        self.sampling  = 1
        self.write_aux = False
        self.expand_2d = False
        self.mpi_split = True
        self.grid      = {}

    def dump(self):
        for attr in sorted(dir(self)):
            if not attr.startswith('_'):
                print "%s = %s" % (attr, getattr(self, attr))

    def get_info(self,pkl_file=None,keys=None,debug=False):
        if pkl_file is None: pkl_file = os.path.join(self.outdir,self._pkl_base+str(0).zfill(4))
        if keys is None: keys = self._basic_info

        filestr = pickle.Unpickler(file(pkl_file))
        data1 = filestr.load()
        data2 = filestr.load()
        data = dict(data1.items() + data2.items())

        for name in keys:
            setattr(self, name, data[name])

        num_frames = len(glob(os.path.join(self.outdir,'*'+self._pkl_base+'*')))
        setattr(self, 'num_frames', num_frames)
        
        if debug:
            self.dump()

        setattr(self,'frame_range',range(0,self.num_frames,self.sampling))
        setattr(self,'q_range',range(self.num_eqn))
        setattr(self,'aux_rage',range(self.num_aux))
        setattr(self,'s_range',range(self.num_dim))

        self.delta = np.asarray(self.delta)

        self.data_zero = data
        self._info_loaded = True

        return

    def read_p(self,frame):
        p_outdir = os.path.join(self.outdir,'_p',self._p_base+str(frame).zfill(4))
        v,n,N = self.ptc_read(p_outdir,num_var=self.num_mp)

        return v

    def build_grid_vecs(self):
        self.lower = self.lower + self.delta/2.0
        self.upper = self.delta*self.num_cells
        
        for k,name in enumerate(self.names):
            setattr(self, name, np.arange(self.lower[0],self.upper[0],self.delta[0]))

        return

    def write_split(self,vec,name,frame=0,irange=None,outdir='./_output'):
        if ((irange is None) or (irange=='all')): irange = range(0,len(vec))

        for n in irange:
            vec[n].tofile(os.path.join(outdir,name+str(n)+'.'+str(frame).zfill(4)))

    def split_solution(self,outdir=None,split_q=True,split_aux=True,update_aux=None,debug=False,ptc_split=None):
        if outdir is None: outdir = self.outdir
        if update_aux is None: update_aux = self.update_aux
        if ptc_split is None: ptc_split=self.mpi_split

        sampling = self.sampling

        if debug: plt.figure()

        for frame in self.frame_range:
            if ptc_split:
                if PETSc.COMM_WORLD.rank==0: print frame
                if split_q:
                    ptcname = self._ptc_base+str(frame).zfill(4)
                    ptcfile = os.path.join(outdir,ptcname)
                    self.ptc_split(outdir=outdir,ptcfile=ptcfile,num_var=self.num_eqn,
                        affix=['sol'+str(frame).zfill(4),'s'],
                        irange=self.s_range)
                if split_aux:
                    auxname = self._aux_base+str(frame).zfill(4)
                    auxfile = os.path.join(outdir,auxname)
                    self.ptc_split(outdir=outdir,ptcfile=auxfile,num_var=self.num_aux,
                        affix=['sol_aux'+str(frame).zfill(4),'n'],
                        poynting=False,irange=self.aux_range)
                    if not update_aux: split_aux=False
            else:
                solution = Solution()
                solution.read(frame,path=outdir,file_format='petsc',read_aux=split_aux)
                q = solution.state.get_q_global()
                if self.num_dim==2:
                    if self.expand_2d:
                        q = q[:,::sampling,::sampling,np.newaxis]
                else:
                    q = q[:,::sampling,::sampling,::sampling]

                if debug:
                    plt.subplot(3,1,1)
                    plt.pcolor(q[1,:,:,16])
                    plt.subplot(3,1,2)
                    plt.pcolor(q[1,:,:,0])

                self.write_split(q,'q',frame=frame,irange=self.q_range,outdir=outdir)

                if split_aux:
                    aux = solution.state.get_aux_global()
                    if self.num_dim==2:
                        if self.expand_2d:
                            aux = aux[:,::sampling,::sampling,np.newaxis]
                    else:
                        aux = aux[:,::sampling,::sampling,::sampling]

                    if debug:
                        plt.subplot(3,1,3)
                        plt.pcolor(aux[0,:,:,16])
                        plt.show()

                    self.write_split(aux,'aux',frame=frame,irange=self.aux_range,outdir=outdir)
                    if not update_aux: split_aux=False

    def ptc_read(self,ptcfile,num_var=1):
        x = PETSc.Vec().create()
        x.setBlockSize(num_var)

        viewer = PETSc.Viewer().createBinary(ptcfile, PETSc.Viewer.Mode.READ)
        x.load(viewer)
        viewer.destroy()

        n, N = x.getSizes()

        v = x.array
        del x
        return v,n,N

    def ptc_split(self,outdir='./',ptcfile='claw.ptc0000',num_var=6,affix=['ptc','A'],debug=True,poynting=True,irange=range(6)):
        x = PETSc.Vec().create()
        x.setBlockSize(num_var)

        viewer = PETSc.Viewer().createBinary(ptcfile, PETSc.Viewer.Mode.READ)
        x.load(viewer)
        viewer.destroy()

        n, N = x.getSizes()

        v = x.array
        v.shape = (-1, num_var)

        z = PETSc.Vec().createMPI((n//num_var, N//num_var))

        if poynting:
            v=np.cross(v[:,:3], v[:,3:])
            num_var = 3

        if (debug and PETSc.COMM_WORLD.rank==0):
            print v.shape
            print ptcfile

        for i in irange:
            mpifile = os.path.join(outdir,affix[1]+str(i).zfill(2)+'.'+affix[0])
            viewer = PETSc.Viewer().createMPIIO(mpifile, PETSc.Viewer.Mode.WRITE)
            z.array = v[:,i]
            z.view(viewer)
            viewer.destroy()

        return

    def set_dir(self,outdir=None):
        if outdir is None:
            outdir = self.savedir
        if not os.path.exists(outdir): os.makedirs(outdir)
        return outdir

class Plot2D(toolbox):
    def __init__(self):
        self.options = {}
        self.savedir = './_plots'
        self.plots   = ['pcolor','quiver','stream']
        super(Plot2D,self).__init__()

    def calc_F(self,u,v):
        F = np.sqrt(u**2 + v**2)
        return F

    def plot_field(self,x,y,u=None,v=None,F=None,contour=False,outdir=None,plot='quiver',figname='_field',format='eps'):
        outdir = self.set_dir(outdir)
        p = 64
        if F is None: F=self.calc_F(u,v)

        plt.close('all')
        plt.figure()
        #fig, axes = plt.subplots(nrows=1)
        if contour:
            plt.hold(True)
            plt.contourf(x,y,F)
        if plot=='quiver':
            plt.quiver(x[::p],y[::p],u[::p],v[::p],scale=0.1)
        
        if plot=='pcolor':
            plt.pcolormesh(x[::4],y[::4],F[::4],cmap=plt.cm.Pastel1)
            plt.colorbar()

        if plot=='stream':
            speed = F[::16]
            plt.streamplot(x[::16], y[::16], u[::16], v[::16], density=(1,1),color='k')

        plt.xlabel('$x$ (a.u.)')
        plt.ylabel('$y$ (a.u.)')
        
        plt.savefig(os.path.join(outdir,figname+'.'+format),format=format,dpi=320,bbox_inches='tight')
        plt.close()

    def make_plots(self,plot_p=True):
        x,y = np.meshgrid(self.x,self.y)
        for frame in self.frame_range:
            print 'ploting frame ',frame
            v = self.read_p(frame)
            v = np.reshape(v,[self.num_mp,self.num_cells[0],self.num_cells[1]], order='F')
            for plot in self.plots:
                print plot
                self.plot_field(x,y,u=v[0].T,v=v[1].T,F=v[3].T,plot=plot,figname=plot+'.'+str(frame).zfill(4),format='png')
            del v

def split_and_xmf(outdir='./_output',update_aux=False,q_range=None,aux_range=None,expand_2d=False):
    from visualization import toolbox
    workbench = toolbox()
    workbench.outdir = outdir
    workbench.read_aux = True
    workbench.update_aux = update_aux
    workbench.get_info()
    workbench.q_range=range(6)
    workbench.aux_range=[0]
    if workbench.num_dim==2:
        workbench.expand_2d=expand_2d
        workbench.q_range=[0,1,2]
        workbench.aux_range=[0,1,2]

    workbench.split_solution()
    if workbench.mpi_split:
        if PETSc.COMM_WORLD.rank==0: workbench.dump()

if __name__=="__main__":
    import sys
    from clawpack.pyclaw import util
    import sys
    args,app_args = util._info_from_argv(sys.argv)
    print app_args
    split_and_xmf(**app_args)
