import os
from glob import glob
import numpy as np
import matplotlib
# matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import pylab as plt
from clawpack.pyclaw import Solution
from scipy.io import loadmat,savemat
from matplotlib.streamplot import  streamplot
matplotlib.rcParams.update({'font.size': 10})

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
    for n in range(5,len(q),sampling):
        if cmap is not None:
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

def poly_plot(verts,cmap=None,sampled=None,view_angles=[10,270],outdir='./',outname='vert_',sampling=10,format='eps'):
    plt.figure()
    ax   = plt.gca(projection='3d')
    if cmap is not None:
        poly = PolyCollection(verts[::sampling],facecolors=cmap[::sampling])
    else:
        poly = PolyCollection(verts[::sampling])

    poly.set_alpha(0.7)
    poly.set_edgecolor('none')
    ax.add_collection3d(poly, zs=sampled[:,0], zdir='y')

    if sampled is not None:
        tticks = np.arange(0,sampled[:,0].max()+1,sampled[:,0].max()/2)
        xticks = np.arange(0,301,100)
        zticks = np.arange(0,sampled[:,6].max()+1,sampled[:,6].max()/2)

        ax.set_yticks(tticks)
        ax.set_xticks(xticks)

        ax.set_xlim3d(0, 300.0)
        ax.set_ylim3d(0, sampled[:,0].max())
        ax.set_zlim3d(0, sampled[:,6].max())

    ax.set_xlabel('$x$',linespacing=3.2)
    ax.set_ylabel('$t$',linespacing=3.2)
    ax.set_zlabel('$I$',linespacing=3.2)
    
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.zaxis.set_ticks_position('none')
    
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    plt.draw()

    for v,view in enumerate(view_angles):
        ax.view_init(elev=view[0], azim=view[1])
        fig_name = outname+str(view[0])+'_'+str(view[1])
        plt.savefig(os.path.join(outdir,fig_name+'.'+format),format=format,dpi=320,bbox_inches='tight')
    return

def assemble_q(path='./_output',frame_plot_range=[0],vecmagnitude=True,poynting=True,poly_verts=True,read_aux=False):
    # create instance of solution object
    num_frames = len(frame_plot_range)
    solution = Solution()
    Q_map_temp = []
    derived_quantities =  {}

    if read_aux: aux_map_temp = []
    if vecmagnitude: I_map_temp = []
    if poynting: S_map_temp = []
    if poly_verts: verts = []
    
    sampled = np.zeros([len(frame_plot_range),7])
    
    # load the frames and assemble values
    for f,frame in enumerate(frame_plot_range):

        solution = Solution()
        solution.read(frame,path=path,file_format='petsc',read_aux=False)

        x = solution.state.grid.x.centers;
        q = solution.state.get_q_global()

        for n,qn in enumerate(q):
            Q_map_temp = np.append(Q_map_temp,qn)

        sampled[f,0] = solution.t
        sampled[f,1] = x[q[0]==q[0].max()][0]
        sampled[f,3] = q[0].max()
        sampled[f,4] = q[1].max()

        if poynting:
            S = q[0]*q[1]
            S_map_temp = np.append(S_map_temp,S,axis=0)

        if vecmagnitude:
            I = np.sqrt(q[0]**2 + q[1]**2)
            sampled[f,6] = I.max()
            I_map_temp = np.append(I_map_temp,I,axis=0)

        if poly_verts:
            verts.append(list(zip(x,I)))

    Q = Q_map_temp.reshape((num_frames,len(q),Q_map_temp.size/(len(q)*num_frames)))

    if vecmagnitude:
        I = (I_map_temp.reshape((num_frames,I_map_temp.size/num_frames)))
        derived_quantities['I'] = I
    
    if poynting:
        S = (S_map_temp.reshape((num_frames,S_map_temp.size/num_frames)))
        derived_quantities['S'] = S

    if poly_verts:
        derived_quantities['vertices'] = verts

    derived_quantities['x'] = x
    derived_quantities['t'] = sampled[:,0]
    derived_quantities['sampled'] = sampled

    return Q,num_frames,derived_quantities

def postprocess_1d(outdir='./_output',multiple=False,overwrite=False,sampling=5,velocity=True,save_mat=True):
    if multiple:
        outdir = outdir+'*'

    outdirs = sorted(glob(outdir))
    print outdirs

    for dirs in outdirs:
        print dirs
        
        figspath = os.path.join(dirs,'_figures')
        binpath  = os.path.join(dirs,'_bin')

        if not os.path.exists(figspath): os.makedirs(figspath)
        if not os.path.exists(binpath): os.makedirs(binpath)
        
        Q,num_frames,derived_quantities = assemble_q(path=dirs,vecmagnitude=True,poynting=True,poly_verts=True,frame_plot_range=frame_plot_range)
            
        colores = get_cmap(num_colors=num_frames)
        
        x = derived_quantities['x']
        sampled = derived_quantities['sampled']

        waterfall_plot(Q[:,0,:],x,sampling=sampling,cmap=colores,num_colors=num_frames,outdir=figspath,outname='waterfall_q0_s'+str(sampling),
                cbar_label='$|q^0|\quad (a.u.)$')
        waterfall_plot(Q[:,1,:],x,sampling=sampling,cmap=colores,num_colors=num_frames,outdir=figspath,outname='waterfall_q1_s'+str(sampling),
                cbar_label='$|q^1|\quad (a.u.)$')

        view_angles = [18,240,0,270,10,270,18,80]
        view_angles = np.reshape(view_angles,(len(view_angles)/2,2))
        poly_plot(verts=derived_quantities['vertices'],sampled=derived_quantities['sampled'],
                cmap=colores,outdir=figspath,view_angles=view_angles,sampling=sampling)

        # calculate velocity
        if velocity:
            sampled = derived_quantities['sampled']
            sampled[1:,2] = (sampled[1:,1]-sampled[:-1,1])/(sampled[1:,0]-sampled[:-1,0])
            sampled = sampled[1::sampling,:]

            plt.figure()
            f, axarr = plt.subplots(2, 1,sharex=True)
            axarr[0].plot(sampled[:,0],sampled[:,3],'o:')
            axarr[0].set_ylabel('$q_{max}\quad (a.u.)$')
            axarr[1].plot(sampled[1:,0],sampled[1:,2],'o:')
            axarr[1].set_ylabel('$v\quad (a.u.)$')
            axarr[1].set_xlabel('time $t\quad (ct)^{-1}$')
            
            plt.draw()
            fig_name = 'velocity'
            plt.savefig(os.path.join(figspath,fig_name+'.eps'),format='eps',dpi=320,bbox_inches='tight')
            plt.close()

        if save_mat:
            summary={'Q':Q,'derived':derived_quantities}
            savemat(os.path.join(figspath,'summary'),summary)

if __name__ == "__main__":
    from clawpack.pyclaw import util
    import sys
    args,app_args = util._info_from_argv(sys.argv)
    print app_args
    postprocess_1d(**app_args)


    

    