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
frame_plot_range = range(1,91,1)

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

def waterfall_plot(q,x,a=None,sampling=10,cmap=None,num_colors=100,outdir='./',outname='waterfall',format='eps',cbar_label='$|q| (a.u.)$'):
    plt.figure()
    fig = plt.figure(facecolor='white')
    ax = plt.axes(frameon=False)
    ax.set_frame_on(False)
    ax.get_xaxis().tick_bottom()
    ax.axes.get_yaxis().set_visible(False)
    ax.hold(True)

    colorVal = '0.3'
    vmax = q[:,:].max()
    for n in range(5,len(q),sampling):
        if cmap is not None:
            colorVal = get_color(value=q[n,:].max(),cmap=cmap,vmax=vmax+.1,num_colors=num_colors)

        ax.plot(x,q[n,:]+n/2.0,label=str(n),color=colorVal,alpha=0.7)
        if a is not None:
            ax.plot(x,10.0*(a[n,:]-1.5)+n/2.0,label=str(n),color='g',alpha=0.7,linewidth=1.0)

    for tic in ax.yaxis.get_major_ticks():
        tic.tick1On = tic.tick2On = False
        tic.label1On = tic.label2On = False

    if cmap is not None:
        scalar = get_smap(vmax=q[:,:].max()+.1,num_colors=sampling)
        cbar = plt.colorbar(scalar)
        cbar.set_label(cbar_label)
    xmin, xmax = ax.get_xaxis().get_view_interval()
    ymin, ymax = ax.get_yaxis().get_view_interval()
    ax.add_artist(Line2D((xmin, xmax), (ymin, ymin), color='black', linewidth=1.5))
    ax.set_xlabel('$x\quad (a.u.)$')
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

def assemble_q(path='./_output',frame_plot_range=[0],vecmagnitude=True,poynting=True,poly_verts=True):
    # create instance of solution object
    num_frames = len(frame_plot_range)
    solution = Solution()
    Q_map_temp = []
    aux_map_temp = []
    derived_quantities = {}
    if vecmagnitude: I_map_temp = []
    if poynting: S_map_temp = []
    if poly_verts: verts = []
    
    sampled = np.zeros([len(frame_plot_range),12])
    
    # load the frames and assemble values
    for f,frame in enumerate(frame_plot_range):

        solution = Solution()
        solution.read(frame,path=path,file_format='petsc',read_aux=True)

        x = solution.state.grid.x.centers;
        q = solution.state.get_q_global()
        aux = solution.state.get_aux_global()[0:2,:]

        for n,qn in enumerate(q):
            Q_map_temp = np.append(Q_map_temp,qn)

        for n,auxn in enumerate(aux):
            aux_map_temp = np.append(aux_map_temp,auxn)

        sampled[f,0] = solution.t
        if len(x[q[0]==q[0].max()])==1:
            it = np.sqrt(q[0]**2 + q[1]**2)
            # print x[np.where(it==it.max())]
            sampled[f,1] = x[np.where(it==it.max())]

        if len(x[q[0]==q[0].max()])>=2:
            sampled[f,1] = x[q[0]==q[0].max()][-1]
            print len(x[q[0]==q[0].max()])
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
    A = aux_map_temp.reshape((num_frames,len(aux),aux_map_temp.size/(len(aux)*num_frames)))

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

    return Q,A,num_frames,derived_quantities

def postprocess_1d(outdir='./_output',base_name='_res_',outsuffix='',read_aux=True,multiple=False,overwrite=False,sampling=5,velocity=True,save_mat=True,poly=False,color=False,lorentz=False):
    if multiple:
        outdir = outdir+'*'+outsuffix

    outdirs = sorted(glob(outdir))
    summary = {}
    print outdirs
    summarypath = '/simdesk/sandbox/emclaw/results/nonlinear/norip/_summary'+outsuffix+'_b_'
    for k,dirs in enumerate(outdirs):
        print dirs
        figspath = os.path.join(dirs,'_figures')
        binpath  = os.path.join(dirs,'_bin')
        ##### alternative for rip/norip test
        base_name_dir = dirs.split('_')
        print base_name_dir
        base_name = base_name_dir[-1]+'_'
        #base_name = base_name_dir[-1]+'_'+base_name_dir[-2]+'_'
        print base_name
        figspath = os.path.join('/simdesk/sandbox/emclaw/results/nonlinear/norip/',base_name_dir[-2])
        print figspath
        #vrip = float(base_name_dir[3].split('v')[1])/100.0
        #if vrip==6349/100.0: vrip =0.6349
        vrip = 1.0/np.sqrt(1.5)
        print vrip
        # base_name = base_name_dir[-1]+'_' #base_name_dir[3]+'_'+base_name_dir[4]+'_'
        # #### alternative for sint
        # base_name_dir = dirs.split('_')
        # print base_name_dir
        # base_name = base_name_dir[-1]+'_'
        # figspath = os.path.join('/simdesk/sandbox/emclaw/results/1D/_sint',base_name_dir[-1])
        # vrip = 0.5
        if not os.path.exists(figspath): os.makedirs(figspath)
        if not os.path.exists(binpath): os.makedirs(binpath)

        Q,A,num_frames,derived_quantities = assemble_q(path=dirs,vecmagnitude=False,poynting=False,poly_verts=False,frame_plot_range=frame_plot_range)

        x = derived_quantities['x']
        sampled = derived_quantities['sampled']

        if color:
            colores = get_cmap(num_colors=num_frames)
        else:
            colores = None

        waterfall_plot(Q[:,0,:],x,sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_q0',
                cbar_label='$|q^0|_{max}\quad (a.u.)$')
        waterfall_plot(Q[:,1,:],x,sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_q1',
                cbar_label='$|q^1|_{max}\quad (a.u.)$')
        waterfall_plot(Q[:,1,:]*Q[:,0,:],x,sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_s',
                cbar_label='$|S|_{max}\quad (a.u.)$')
        waterfall_plot(np.sqrt(Q[:,1,:]**2 + Q[:,0,:]**2),x,sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_i',
                cbar_label='$I_{max}\quad (a.u.)$')
        if read_aux:
            print 'read'
            waterfall_plot(10.0*np.sqrt(A[:,1,:]*A[:,0,:]),x,sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_n',
                    cbar_label='$n_{max}\quad (a.u.)$')
            waterfall_plot(Q[:,1,:]*Q[:,0,:],x,np.sqrt(A[:,1,:]*A[:,0,:]),sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_sn',
                    cbar_label='$|S|_{max}\quad (a.u.)$')
            waterfall_plot(np.sqrt(Q[:,1,:]**2 + Q[:,0,:]**2),x,np.sqrt(A[:,1,:]*A[:,0,:]),sampling=5,cmap=colores,num_colors=num_frames,outdir=figspath,outname=base_name+'waterfall_in',
                    cbar_label='$I_{max}\quad (a.u.)$')
        if poly:
            view_angles = [18,240,0,270,10,270,18,80]
            view_angles = np.reshape(view_angles,(len(view_angles)/2,2))
            poly_plot(verts=derived_quantities['vertices'],sampled=derived_quantities['sampled'],
                    cmap=colores,outdir=figspath,view_angles=view_angles,sampling=5,outname=base_name+'vert_')

        # calculate velocity
        if velocity:
            sampled = derived_quantities['sampled']
            sampled[:,2] = np.gradient(sampled[:,1],sampled[2,0]-sampled[1,0])
            sampled[:,7] = np.sqrt(sampled[:,3]**2 + sampled[:,4]**2)
            sampled[:,8] = np.gradient(sampled[:,7],sampled[2,0]-sampled[1,0])
            sampled[:,9] = sampled[:,3]*sampled[:,4]
            sampled[:,10] = np.gradient(sampled[:,9],sampled[2,0]-sampled[1,0])

            tt = sampled[6::sampling,0]
            xx = sampled[6::sampling,1]
            dt = tt[1] - tt[0]
            s  = sampled[6::sampling,9]
            i  = sampled[6::sampling,7]
            dx = np.gradient(sampled[6::sampling,1],dt)
            di = np.gradient(sampled[6::sampling,7],dt)
            ds = np.gradient(sampled[6::sampling,9],dt)
            bis = [tt,xx,dx,di,ds,s,i]
            summary['outdir'+str(k)] = dirs
            summary['vrip'+str(k)]    = vrip
            summary['sampled'+str(k)] = sampled
            summary['dt'+str(k)]  = dt
            summary['bis'+str(k)] = bis

            if lorentz:
                gamma = 1.0/np.sqrt(1 - vrip**2)
                xp  = gamma*((xx) - vrip*tt)
                tp  = gamma*(tt - vrip*(xx))
                vp = (xp[1:] - xp[0:-1])/(tp[1:] - tp[0:-1])
                dxp = np.gradient(xp,tp)

                sampled_lorentz = [tp,xp,vp,dxp]
                summary['lorentz'+str(k)] = sampled_lorentz
                
                labels = ['$t^p\quad (a.u.)$','$x^p_{max}\quad (a.u.)$','$dx^p_{max}/dt^p\quad (a.u.)$','$dx^p_{max}/dt^p\quad (a.u.)$']
                savens = ['x_bis','v_bis','xp_bis']
                for i in range(1,len(sampled_lorentz)):
                    tpl = tp
                    if i==2:
                        tpl = tp[1:]
                    plot_single(tpl,sampled_lorentz[i],xlabel=labels[0],ylabel=labels[i],
                        figspath=figspath,figname=base_name+'lorentz_'+savens[i-1])

            plot_single(sampled[6::sampling,0],sampled[6::sampling,1],ylabel='$x_{max}\quad (a.u.)$',
                figspath=figspath,figname=base_name+'x')

            plot_single(sampled[6::sampling,0],sampled[6::sampling,7],ylabel='$I_{max}\quad (a.u.)$',
                figspath=figspath,figname=base_name+'i')

            plot_single(sampled[6::sampling,0],sampled[6::sampling,9],ylabel='$|S|_{max}\quad (a.u.)$',
                figspath=figspath,figname=base_name+'s')

            plot_single(tt,dx,ylabel='$dx_{max}/dt\quad (a.u.)$',
                figspath=figspath,figname=base_name+'dxdt_bis',ylim=[0.4,0.8],vrip=vrip)

            plot_single(tt,xx,ylabel='$x_{max}\quad (a.u.)$',
                figspath=figspath,figname=base_name+'x_bis',xrip=vrip*tt+25.0)

            plot_single(tt,di,ylabel='$dI_{max}/dt\quad (a.u.)$',
                figspath=figspath,figname=base_name+'didt_bis')

            plot_single(tt,ds,ylabel='$d|S|_{max}/dt\quad (a.u.)$',
                figspath=figspath,figname=base_name+'dsdt_bis')

    if not os.path.exists(summarypath): os.makedirs(summarypath)
    ndirs = len(outdirs)
    plot_summary(0,1,summary,ndirs=ndirs,ylabel='$x_{max}\quad (a.u.)$',
                figspath=summarypath,figname='x_bis')
    plot_summary(0,2,summary,ndirs=ndirs,ylabel='$dx_{max}/dt\quad (a.u.)$',
                figspath=summarypath,figname='dxdt_bis',ylim=[0.4,0.8])
    plot_summary(0,3,summary,ndirs=ndirs,ylabel='$dI_{max}/dt\quad (a.u.)$',
                figspath=summarypath,figname='didt_bis')
    plot_summary(0,4,summary,ndirs=ndirs,ylabel='$d|S|_{max}/dt\quad (a.u.)$',
                figspath=summarypath,figname='dsdt_bis')
    plot_summary(0,5,summary,ndirs=ndirs,ylabel='$|S|_{max}\quad (a.u.)$',
                figspath=summarypath,figname='s_bis')
    plot_summary(0,6,summary,ndirs=ndirs,ylabel='$I_{max}\quad (a.u.)$',
                figspath=summarypath,figname='i_bis')
    if lorentz:
        plot_summary(0,1,summary,ndirs=ndirs,dictsrc='lorentz',
            xlabel='$t^p\quad (a.u.)$',ylabel='$x^p_{max}\quad (a.u.)$',
            figspath=summarypath,figname='lorentz_x_bis')
        plot_summary(0,2,summary,p=1,ndirs=ndirs,dictsrc='lorentz',
            xlabel='$t^p\quad (a.u.)$',ylabel='$dx^p_{max}/dt^p\quad (a.u.)$',
            figspath=summarypath,figname='lorentz_vp_bis')
        plot_summary(0,3,summary,ndirs=ndirs,dictsrc='lorentz',
            xlabel='$t^p\quad (a.u.)$',ylabel='$dx^p_{max}/dt^p\quad (a.u.)$',
            figspath=summarypath,figname='lorentz_dxdt_bis')

    if save_mat:
        savemat(os.path.join(figspath,base_name+'summary'),summary)

def plot_summary(x,y,dictionary,ndirs=4,p=0,dictsrc='bis',label_base='$v^{0}=$',xlabel='$t\quad (ct)^{-1}$',ylabel='y',shape='--',figspath='./_output',figname='figure',ylim=None):
    plt.close('all')
    plt.figure()
    f, ax = plt.subplots(1, 1,sharex=True)
    ax.hold(True)
    chival=[0.55,0.59,0.60,0.61,0.63]
    #chival=[4,12,24,128]
    #chival=[0.1,0.01,0.001,0.0]
    for k in range(ndirs):
        xx = dictionary[dictsrc+str(k)][x][p:]
        yy = dictionary[dictsrc+str(k)][y]
        ax.plot(xx,yy,shape,label=label_base+' $'+str(chival[k])+'$')

    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    if ylim is not None:
        ax.set_ylim(ylim)
    ax.legend(frameon=False,loc='upper center', bbox_to_anchor=(0.5, -0.15),ncol=len(chival))
    plt.draw()
    plt.savefig(os.path.join(figspath,figname+'.eps'),format='eps',dpi=320,bbox_inches='tight')
    plt.close()

def plot_together(x,y1,y2,xlabel='$t\quad (ct)^{-1}$',y1label='y1',y2label='y2',shape='.:',figspath='./_output',figname='figure'):
        plt.close('all')
        plt.figure()
        f, axarr = plt.subplots(2, 1,sharex=True)
        axarr[0].plot(x,y1,shape)
        axarr[0].set_ylabel(y1label)
        axarr[1].plot(x,y2,shape)
        axarr[1].set_ylabel(y2label)
        axarr[1].set_xlabel(xlabel)

        plt.draw()
        plt.savefig(os.path.join(figspath,figname+'.eps'),format='eps',dpi=320,bbox_inches='tight')
        plt.close()

def plot_single(x,y,xlabel='$t\quad (ct)^{-1}$',ylabel='y',shape='--',figspath='./_output',figname='figure',ylim=None,vrip=None,xrip=None):
        plt.close('all')
        plt.figure()
        f, axarr = plt.subplots(1, 1,sharex=True)
        axarr.plot(x,y,shape)
        if vrip is not None:
            axarr.plot(x,vrip*np.ones(len(x)),'r')
        if xrip is not None:
            axarr.plot(x,xrip,'r')
        axarr.set_ylabel(ylabel)
        axarr.set_xlabel(xlabel)
        if ylim is not None:
            axarr.set_ylim(ylim)
        plt.draw()
        plt.savefig(os.path.join(figspath,figname+'.eps'),format='eps',dpi=320,bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    from clawpack.pyclaw import util
    import sys
    args,app_args = util._info_from_argv(sys.argv)
    print app_args
    postprocess_1d(**app_args)
