fontsize = 18
import os
from glob import glob
import numpy as np
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

# plot frame_plot_range together, define range
frame_plot_range = range(1,101,1)

def get_cmap(colormap='jet',num_colors=100):
    cmap = cm  = plt.get_cmap(colormap)
    cNorm     = colors.Normalize(vmin=0, vmax=num_colors)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    colores   = [scalarMap.to_rgba(k) for k in range(0,num_colors) ]    
    return colores

def get_smap(colormap='jet',vmin=0.0,vmax=10.0,num_colors=100):
    cmap = cm  = plt.get_cmap(colormap)
    cNorm     = colors.Normalize(vmin=0, vmax=vmax)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    scalarMap._A = []
    return scalarMap

def get_color(value,cmap,vmin=0.0,vmax=10.0,num_colors=100):
    values    = np.arange(0,vmax+vmax/num_colors,vmax/num_colors)
    colorVal  = cmap[np.where(values>=value)[0][0]]
    return colorVal

def waterfall_plot(q,x,sampling=10,cmap=None,num_colors=100,outdir='./',outname='waterfall',format='eps',cbar_label='$|q| (a.u.)$'):
    plt.figure()
    plt.hold(True)
    colorVal = 'b'
    vmax = q[:,:].max()
    # print vmax,len(q)
    for n in range(0,len(q),sampling):
        if cmap is not None:
            # print q[n,:].max()
            colorVal = get_color(value=q[n,:].max(),cmap=cmap,vmax=vmax+.1,num_colors=num_colors)

        plt.plot(x,q[n,:]+n/10.0,label=str(n),color=colorVal,alpha=0.7)
    ax = plt.gca()
    for tic in ax.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False

    if cmap is not None:
        scalar = get_smap(vmax=q[:,:].max()+.1,num_colors=sampling)
        cbar = plt.colorbar(scalar)

    plt.xlabel('$x\quad (a.u.)$')
    cbar.set_label(cbar_label)
    plt.draw()

    plt.savefig(os.path.join(outdir,outname+'.'+format),format=format,dpi=320,bbox_inches='tight')
    plt.close()
    return

def get_grid(solution,sampling):
    x = solution.state.grid.c_centers[0][::sampling[0],::sampling[1]]
    y = solution.state.grid.c_centers[1][::sampling[0],::sampling[1]]
    x_lower,x_upper = solution.state.grid.dimensions[0].lower,solution.state.grid.dimensions[0].upper
    y_lower,y_upper = solution.state.grid.dimensions[1].lower,solution.state.grid.dimensions[1].upper
    mx,my = solution.state.grid.dimensions[0].num_cells,solution.state.grid.dimensions[1].num_cells
    dimension = np.asarray([[x_lower,y_lower],[x_upper,y_upper],[mx,my]])
    return x,y,dimension

def read_q(path='./_output',frame=0,read_aux=False,read_grid=False,read_q=True,qi=[0,1,2],sampling=[1,1],loc_max=True):
    solution = Solution()
    
    solution.read(frame,path=path,file_format='petsc',read_aux=read_aux)
    sol_compilation = {}
    sol_compilation['sampling'] = sampling

    if read_grid:
        x,y,dimension = get_grid(solution,sampling)
        dict2 = {'x':x,'y':y,'nodes':dimension[0:2,:],'num_cells':dimension[-1,:]}
        sol_compilation.update(dict2)
    
    if read_aux:
        sol_compilation.update({'aux':solution.state.get_aux_global()[0:3,::sampling[0],::sampling[1]]})

    if read_q:
        q_all = solution.state.get_q_global()[:,::sampling[0],::sampling[1]]
        for qn in qi:
            q = solution.state.get_q_global()[qn,::sampling[0],::sampling[1]]
            sol_compilation.update({'q'+str(qn):q})
            if loc_max:
                if not read_grid:
                    x,y,dimension = get_grid(solution,sampling)
                indmax = q[:,:]==q[:,:].max()
                sol_compilation.update({'indmax'+str(qn):indmax})
                sol_compilation.update({'xmax':x[indmax],'ymax':y[indmax]})

    sol_compilation.update({'t':solution.t})        

    return sol_compilation if not read_q else q_all,sol_compilation

def extract_cut_line(q,x,y,plot=False,dim='xy',points=[0.0,5.0,199.0,5.0],num_points=1000,outdir='./',outname='cut',format='eps'):
    if not os.path.exists(outdir): os.makedirs(outdir)
    import scipy.ndimage

    # get nodes based on [xi,yi],[xf,yf]
    points = np.asarray(points)
    points = points.reshape((2,2))

    # sanity check
    if x[:,0].max()<=points[1,0]:
        points[1,0] = x[:,0].max() - 2.0*(x[1,0]-x[0,0])

    if y[0,:].max()<=points[0,1]:
        points[0,1] = y[0,:].max() - 2.0*(y[0,1]-y[0,0])

    if y[0,:].max()<=points[1,1]:
        points[1,1] = y[0,:].max() - 2.0*(y[0,1]-y[0,0])

    xnodes = [np.where(x[:,0]>=(points[i,0]))[0][0] for i in range(0,2)]
    
    ynodes = [np.where(y[0,:]>=(points[i,1]))[0][0] for i in range(0,2)]

    temp_flag = True
    
    if (ynodes[0]-ynodes[1])==0.0:
        qi = q[:,ynodes[0]]
        di = x[:,ynodes[0]]
        temp_flag = False

    if (xnodes[0]-xnodes[1])==0.0:
        qi = q[xnodes[0],:]
        di = y[xnodes[0],:]
        temp_flag = False

    if temp_flag:
        # create an array with distances
        d = np.sqrt(x**2+y)

        nodes = np.append(xnodes,ynodes,axis=0)
        nodes = nodes.reshape(2,2).T

        xn,yn = np.linspace(nodes[0,0], nodes[1,0], num_points), np.linspace(nodes[0,1], nodes[1,1], num_points) 

        # Extract the values along the line, using cubic interpolation
        qi = scipy.ndimage.map_coordinates(q, np.vstack((xn,yn)))
        di = scipy.ndimage.map_coordinates(d, np.vstack((xn,yn)))

    if plot:
        plt.close('all')
        plt.figure()
        fig, axes = plt.subplots(nrows=2)
        axes[0].pcolormesh(x,y,q,cmap='jet',vmin=0,vmax=1.2)
        axes[0].plot(points[:,0],points[:,1], 'r.:')
        axes[0].set_xlabel('$x$')
        axes[0].set_ylabel('$y$')

        axes[1].plot(di,qi)
        axes[1].set_xlabel('$d$')
        axes[1].set_ylabel('$q$')

        plt.savefig(os.path.join(outdir,outname+'.'+format),format=format,dpi=320,bbox_inches='tight')
        plt.close()
    
    return di,qi

def plot_field(x,y,u,v,contour=True,outdir='./',outname='field',format='eps',cut_long=False):
    if not os.path.exists(outdir): os.makedirs(outdir)
    if cut_long: nrows=2
    plt.close('all')
    plt.figure()
    fig, axes = plt.subplots(nrows=nrows)
    if contour:
        axes[0].hold(True)
        axes[0].contour(x,y,np.sqrt(u**2 + v**2))
    axes[0].quiver(x,y,u,v,angles='xy', scale_units='xy', scale=10)
    axes[0].set_xlabel('$x$')
    axes[0].set_ylabel('$y$')
    if cut_long:
        axes[1].hold(True)
        axes[1].plot(x[:,x.shape[1]/2],u[:,x.shape[1]/2],':',label='Fx',color='r')
        axes[1].plot(x[:,x.shape[1]/2],v[:,x.shape[1]/2],':',label='Fy',color='g')
        axes[1].set_xlabel('$x$')
        axes[1].set_ylabel('$|field|$')

    plt.savefig(os.path.join(outdir,outname+'.'+format),format=format,dpi=320,bbox_inches='tight')
    plt.close()
    return

def plot_streamline(x,y,u,v,outdir='./',outname='field_stream',format='eps'):
    if not os.path.exists(outdir): os.makedirs(outdir)
    speed = np.sqrt(u**2+v**2)
    plt.close('all')
    plt.figure()
    plt.streamplot(x, y, u, v, color=speed,linewidth=5*speed/speed.max(),cmap=plt.cm.jet)
    plt.savefig(os.path.join(outdir,outname+'.'+format),format=format,dpi=320,bbox_inches='tight')
    plt.close()
    return

def Poynting(q,plot=False,streamline=False,x=None,y=None,outdir='./',suffix='',format='eps',cut_long=False):
    u =      q[1]*q[2]
    v = -1.0*q[0]*q[2]

    if plot:
        outname = 'S'+suffix
        plot_field(x,y,u,v,outdir=outdir,outname=outname,format=format,cut_long=cut_long)
    
    if streamline:
        outname = 'stream'+suffix
        plot_streamline(y,x,u,v,outdir=outdir,outname=outname,format=format)    

    return u,v

def assemble_q(path='./_output',frame_plot_range=[0],poynting=True,read_aux=False,qi=[0,1,2],update_aux=False,sampling=[1,1],
    frame_split=False,split_q=False,figspath='./',binpath='./',cut=True):
    # create instance of solution object
    num_frames = len(frame_plot_range)
    solution = Solution()
    Q_map_temp = [[],[]]
    d_map_temp = [[],[]]
    sol = {}
    
    # print frame_plot_range
    sampled = np.zeros([num_frames,9])
    plot_aux = True
    # load the frames and assemble values
    for f,frame in enumerate(frame_plot_range):
        # print frame
        ql = [[],[]]
        dl = [[],[]]
        if f==0:
            q,grid = read_q(path=path,frame=frame,read_q=False,read_grid=True,sampling=sampling)
            sol['grid'] = grid

        if read_aux:
            if f==0 or update_aux:
                q,aux = read_q(path=path,frame=frame,read_q=False,read_aux=True,qi=[1,2],sampling=sampling)
                sol['aux'+str(f)] = aux

                if plot_aux:
                    plt.close('all')
                    plt.figure()
                    fig, axes = plt.subplots(nrows=1)
                    axes.pcolormesh(grid['x'],grid['y'],aux['aux'][1,:,:],cmap='jet')
                    axes.set_xlabel('$x$')
                    axes.set_ylabel('$y$')
                    aux_name = 'aux'+str(frame).zfill(4)
                    plt.savefig(os.path.join(figspath,aux_name+'.png'),format='png',dpi=320,bbox_inches='tight')
                    plt.close()
                    plot_aux = update_aux

        q, solution = read_q(path=path,frame=frame,qi=qi,sampling=sampling)

        for k in xrange(1,3,1):        
            dl[k-1],ql[k-1] = extract_cut_line(q[k],grid['x'],grid['y'],plot=cut,num_points=1000,
                outdir=os.path.join(figspath,'cut_x'),outname='cut_q'+str(k-1)+'_'+str(frame).zfill(4),format='png')


            Q_map_temp[k-1] = np.append(Q_map_temp[k-1],ql[k-1])
            d_map_temp[k-1] = np.append(d_map_temp[k-1],dl[k-1])

        sol['q'+str(f)] = q

        sampled[f,0] = solution['t']
        sampled[f,1] = solution['xmax'].max()
        sampled[f,2] = solution['ymax'].max()
        sampled[f,3] = solution['q1'].max()
        sampled[f,4] = solution['q2'].max()

        if poynting:
            p1 = q.shape[1]/100
            p2 = q.shape[2]/50
            Sx,Sy = Poynting(q[:,::p1,::p2],plot=True,streamline=True,x=grid['x'][::p1,::p2],y=grid['y'][::p1,::p2],
                    outdir=os.path.join(figspath,'Poyinting'),suffix=str(frame).zfill(4),format='png',cut_long=True)
        

        if split_q:
            for qn in range(0,len(q)):
                q_temp = q[qn,:,:,np.newaxis]
                aux_temp = aux['aux'][qn,:,:,np.newaxis]
                q_temp.tofile(os.path.join(binpath,'3D_q'+str(qn)+'.'+str(frame).zfill(4)))
                aux_temp.tofile(os.path.join(binpath,'3D_aux'+str(qn)+'.'+str(frame).zfill(4)))


        if frame_split: 
            sol['Sx'+str(f)] = Sx
            sol['Sy'+str(f)] = Sy
            sol['sampled'+str(f)] = sampled
            savemat(os.path.join(figspath,'sol'+str(frame).zfill(4)),sol)

    Q = np.reshape(Q_map_temp,(2,num_frames,np.shape(Q_map_temp)[1]/num_frames))
    d = np.reshape(d_map_temp,(2,num_frames,np.shape(d_map_temp)[1]/num_frames))


    sol['sampled'] = sampled
    sol['num_frames'] = num_frames

    return Q,d,sol

def postprocess(outdir='./_output',multiple=False,overwrite=False,sampling=5,save_mat=True,
    frame_split=False,split_q=False,update_aux=True,poynting=True,cut=True):
    if multiple:
        outdir = outdir+'*'

    outdirs = sorted(glob(outdir))
    # print outdirs

    for dirs in outdirs:
        # print dirs
        
        figspath = os.path.join(dirs,'_figures')
        binpath  = os.path.join(dirs,'_bin')

        temp_flag = True
        if not os.path.exists(figspath):
            os.makedirs(figspath)
        else:
            temp_flag = False
        
        if not os.path.exists(binpath):
            os.makedirs(binpath)
        else:
            temp_flag = False

        if not temp_flag:
            if overwrite: temp_flag = True

        if temp_flag:
            Q,d,summary = assemble_q(path=dirs,frame_plot_range=frame_plot_range,poynting=poynting,read_aux=True,
                update_aux=update_aux,frame_split=frame_split,split_q=split_q,figspath=figspath,binpath=binpath,cut=cut)

            # print Q.shape

            colores = get_cmap(num_colors=summary['num_frames']+1)
            waterfall_plot(Q[0,:,:],d[0,0,:],sampling=sampling,cmap=colores,num_colors=summary['num_frames'],
                outdir=figspath,outname='waterfall_q0_s'+str(sampling),cbar_label='$|q^0|\quad (a.u.)$')
            waterfall_plot(Q[1,:,:],d[1,0,:],sampling=sampling,cmap=colores,num_colors=summary['num_frames'],
                outdir=figspath,outname='waterfall_q1_s'+str(sampling),cbar_label='$|q^1|\quad (a.u.)$')
            
            if save_mat:
                savemat(os.path.join(figspath,'summary'),summary)

# num_frames = len(frame_plot_range)
if __name__ == "__main__":
    # print 'la'
    from clawpack.pyclaw import util
    import sys
    args,app_args = util._info_from_argv(sys.argv)
    # print app_args
    # print frame_plot_range
    postprocess(**app_args)