import os
from glob import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import pylab as plt
from clawpack.pyclaw import Solution
from scipy.io import loadmat,savemat

# plot frame_plot_range together, define range
frame_plot_range = range(1,101,1)

# whether to generate a summary matlab file
genmatfile = True

# pyclaw source path and filenames
testpath = '/simdesk/sandbox/emclaw/maxwell_vc_1d/'
testbase = '_output_'
testnums = range(13,14,1)

# output peraphernalia
figspath = '/simdesk/sandbox/emclaw/analysis/results/postfigs/100pts'

def get_cmap(colormap='jet',num_colors=100):
    cmap = cm  = plt.get_cmap(colormap)
    cNorm     = colors.Normalize(vmin=0, vmax=num_colors)
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cmap)
    colores   = [scalarMap.to_rgba(k) for k in range(0,num_colors) ]    
    return colores

def get_color(value,cmap,vmin=0.0,vmax=10.0,num_colors=100):
    values    = np.arange(0,vmax+vmax/num_colors,vmax/num_colors)
    colorVal  = colores[np.where(values>=value)[0][0]]
    return colorVal

# create folder to save results
if not os.path.exists(figspath): os.makedirs(figspath)

# superposition of solution at given frame_plot_range
sampled = np.zeros([len(testnums),len(frame_plot_range),7])

# plot peraphernalia
matplotlib.rcParams.update({'font.size': 10})
colores = get_cmap()

# matlab file for I_map
I_map = {}
Q_map = {}

for k,enddir in enumerate(testnums):
    # get the directory
    dirs = os.path.join(testpath,testbase+str(enddir))
    
    # create instance of solution object
    solution = Solution()

    I_map_temp = []
    verts = []
    colorVal = []
    time = []
    # load the frames
    for f,frame in enumerate(frame_plot_range):
        solution = Solution()
        solution.read(frame,path=dirs,file_format='petsc',read_aux=False)

        x = solution.state.grid.x.centers;
        q = solution.state.get_q_global()
        
        time = np.append(time,solution.t)
        sampled[k,f,0] = solution.t
        sampled[k,f,1] = x[q[0]==q[0].max()]
        sampled[k,f,3] = q[0].max()
        sampled[k,f,4] = q[1].max()

        colorVal = scalarMap.to_rgba(f)
        
        for n,qn in enumerate(q):
            plt.figure(n)
            plt.hold(True)
            plt.plot(x,qn,':',label=str(frame),color=colorVal)
            plt.draw()

        I = np.sqrt(q[0]**2 + q[1]**2)
        sampled[k,f,6] = I.max()

        I_map_temp = np.append(I_map_temp,I,axis=0)
        verts.append(list(zip(x,I)))

        plt.figure(int(len(q)+1))
        plt.hold(True)
        plt.plot(x,I,':')


    plt.figure(int(len(q)+2))
    ax   = plt.gca(projection='3d')
    poly = PolyCollection(verts,facecolors=colors)
    poly.set_alpha(0.7)
    poly.set_edgecolor('none')
    ax.add_collection3d(poly, zs=sampled[k,:,0], zdir='y')

    tticks = np.arange(0,sampled[k,:,0].max()+1,sampled[k,:,0].max()/2)
    xticks = np.arange(0,301,100)
    zticks = np.arange(0,sampled[k,:,6].max()+1,sampled[k,:,6].max()/2)

    ax.set_yticks(tticks)
    ax.set_xticks(xticks)

    ax.set_xlabel('$x$',linespacing=3.2)
    ax.set_xlim3d(0, 300.0)
    ax.set_ylabel('$t$',linespacing=3.2)
    ax.set_ylim3d(0, sampled[k,:,0].max())
    ax.set_zlabel('$I$',linespacing=3.2)
    ax.set_zlim3d(0, sampled[k,:,6].max())
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.zaxis.set_ticks_position('none')
    
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.draw()

    X,T = np.meshgrid(x,sampled[k,:,0])
    I_map_temp = I_map_temp.reshape((len(frame_plot_range),q.shape[1]))

    # plt.figure(int(len(q)+3))
    # plt.pcolormesh(X,T,I_map_temp,cmap=jet)
    # plt.xlabel('$x (ct)$')
    # plt.ylabel('$t (ct)^{-1}$')
    # plt.colorbar()
    # plt.draw()

    I_map['x'+str(k)] = X
    I_map['t'+str(k)] = T
    I_map['I'+str(k)] = I_map_temp

    for n in range(0,len(q)):
        plt.figure(n)
        plt.axis([0,300,0,5])
        # plt.title('Evloution of $q^'+str(n)+'_n$ at '+str(int(q.shape[1])))
        plt.xlabel('$x\quad (ct)$')
        plt.ylabel('$q\quad (a.u.)$')
        plt.draw()

        fig_name = 'q_'+str(n)+'_'+str(int(q.shape[1]))
        plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    
        plt.close(n)

    plt.figure(int(len(q)+1))
    plt.axis([0,300,0,7])
    # plt.title('Evloution of $I$ at '+str(int(q.shape[1])))
    plt.xlabel('$x\quad (ct)$')
    plt.ylabel('$I\quad (a.u.)$')
    plt.draw()

    fig_name = 'I_'+str(int(q.shape[1]))
    plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    plt.close(int(len(q)+1))
    
    plt.figure(int(len(q)+2))
    view_angles = [18,240,0,270,10,270,18,80]
    view_angles = np.reshape(view_angles,(len(view_angles)/2,2))
    for v,view in enumerate(view_angles):
        ax.view_init(elev=view[0], azim=view[1])
        fig_name = 'vert_'+str(view[0])+'_'+str(view[1])+'_'+str(int(q.shape[1]))
        plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    
    plt.close(int(len(q)+2))

    plt.figure(int(len(q)+3))
    fig_name = 'I_map_ '+str(int(q.shape[1]))
    plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    plt.close(int(len(q)+3))

    # calculate velocity
    sampled[k,1:,2] = (sampled[k,1:,1]-sampled[k,:-1,1])/(sampled[k,1:,0]-sampled[k,:-1,0])
    plt.figure(int(len(q)+4))
    f, axarr = plt.subplots(2, 1,sharex=True)
    axarr[0].plot(sampled[k,:,0],sampled[k,:,3],'o:')
    # axarr[0].set_title('$q_{max}$ vs $t$')
    axarr[0].set_ylabel('$q_{max}\quad (a.u.)$')
    axarr[1].plot(sampled[k,1:,0],sampled[k,1:,2],'o:')
    # axarr[1].set_title('$v$ vs $t$')
    axarr[1].set_ylabel('$v\quad (a.u.)$')
    axarr[1].set_xlabel('time $t\quad (ct)^{-1}$')
    
    plt.draw()
    fig_name = 'v_'+str(int(q.shape[1]))
    plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    plt.close(int(len(q)+4))

savemat(os.path.join(figspath,'I_summary'),I_map)