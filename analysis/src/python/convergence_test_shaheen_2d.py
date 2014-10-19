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
matplotlib.rcParams['lines.markersize'] = 10
matplotlib.rcParams['lines.color'] = 'r'
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib import pylab as plt
from clawpack.pyclaw import Solution
from scipy.io import loadmat,savemat
from matplotlib.streamplot import  streamplot
from scipy.io import loadmat,savemat

# pyclaw source path and filenames
testpath = '/media/shaheen/dev/results/maxwell_2d/convergence_averaged'
testbase = '_output_plane_'

# output peraphernalia
savepath = '/simdesk/sandbox/emclaw/results/2D/convergence_averaged/summary'
figspath = '/simdesk/sandbox/emclaw/results/2D/convergence_averaged/summary/figures'
repopath = '/simdesk/sandbox/emclaw/results/2D/convergence_averaged/summary/reports'
outprefix = 'matgrid_'

# define frame of interest for analytic vs numeric convergence
frame = 45
qn = 2

# local range configuration
local_range_delta = 2.0

# debug flag
debug_flag = True

def debug_convergence_message(results,restype):
    print '\n'+restype+' (absolute)\t\t'+str(results[0])
    print restype+' (relative)\t\t'+str(results[1])

def self_convergence(claw_old,claw_new,dim='xy',debug=False):
    pfinal1 = claw_old.state.get_q_global()
    pfinal2 = claw_new.state.get_q_global()
    grid = claw_new.state.grid

    # average claw new
    if claw_new.state.num_dim==1:
        pfinal2 = (pfinal2[:,::2] + pfinal2[:,1::2])/2.0
    elif claw_new.state.num_dim==2:
        pfinal2 = (pfinal2[:,::2,::2] + pfinal2[:,1::2,::2] + pfinal2[:,::2,:1:2] + pfinal2[:,1::2,1::2])/4.0
    elif claw_new.state.num_dim==3:
        pfinal2 = (pfinal2[:,::2,::2,::2] + pfinal2[:,1::2,::2,::2] + pfinal2[:,::2,:1:2,::2] + pfinal2[:,1::2,1::2,::2] + 
            pfinal2[:,::2,::2,1::2] + pfinal2[:,1::2,::2,1::2] + pfinal2[:,::2,:1:2,1::2] + pfinal2[:,1::2,1::2,1::2])/8.0

    pfinal2 = pfinal2[qn]
    pfinal1 = pfinal1[qn]

    pfinal2 = pfinal2.reshape(-1)
    pfinal1 = pfinal1.reshape(-1)

    self_difference_absolute = np.prod(grid.delta)*np.linalg.norm(pfinal1-pfinal2,ord=1)
    self_difference_relative = self_difference_absolute/np.linalg.norm(pfinal2,ord=1)

    self_difference = [self_difference_absolute, self_difference_relative]

    if debug:
        debug_convergence_message(self_difference,'self')
    
    return self_difference

def reduce_dim(q):
    if len(q.shape)==2:
        q = q[:,q.shape[1]/2]
    elif len(q.shape)==3:
        q = q[:,q.shape[1]/2,q.shape[2]/2]

    return q

def create_dir(directory):
    if not os.path.exists(directory): os.makedirs(directory)
    return

def plot_convergence(x,y,rate=None,log=True,outdir='./',figname='_plot',case='',title='',format='eps',tolatex=False):
    create_dir(outdir)
    log_name  = ''
    log_label = ''
    log_title = ''
    if not tolatex:
        log_title = 'Log Log '
    if log:
        log_name  = 'loglog_'

    plt.figure()
    if rate is not None:
        m = rate[0]
        c = rate[1]
        plt.hold(True)
        plt.loglog(x,10.0**(m*np.log10(x)+c),'r*-.')
    if log:
        plt.loglog(x,y,'bo--')
    else:
        plt.plot(x,y,'o:')
    if not tolatex:
        plt.title(title+log_title+case+' convergence')
    
    plt.xlabel('$'+log_label+'mx$')
    plt.ylabel('$'+log_label+'\Delta q$')

    plt.grid(True)
    plt.draw()
    case = case.replace(" ", "_")
    fig_save_name = figname+'_'+log_name+case+'.'+format
    figpath = os.path.join(outdir,fig_save_name)
    plt.savefig(figpath,format=format,dpi=320,bbox_inches='tight')
    plt.close()

    if tolatex:
        caption = ''
        if log:
            caption = 'Log Log '
        caption = caption + 'plot of '+case.replace("_"," ")+' convergence test'
        caption.capitalize()
        gen_latex_fig(figpath,caption=caption)
    return

def gen_latex_fig(figpath,caption=''):
    create_dir(repopath)
    dstfile = os.path.join(repopath,'summary_figures.tex')
    f = open(dstfile,'a')
    strt = r'\begin{figure}[h!]'
    strt = strt + '\n \t' + r'\centering'
    strt = strt + '\n \t' + r'\includegraphics[width=0.75\textwidth]{'+os.path.relpath(figpath,repopath)+'}'
    strt = strt + '\n \t' + r'\caption{'+caption+'}'
    strt = strt + '\n \t' + r'\label{fig:xxxx}'
    strt = strt + '\n' + r'\end{figure}'
    strt = strt + '\n\n'
    f.write(strt)
    f.close()
    return

def linear_fit(x,y):
    A = np.vstack([x,np.ones(x.size)]).T
    m,c = np.linalg.lstsq(A,y)[0]
    return m,c

# get pyclaw directories for results
testdirs = sorted(glob(os.path.join(testpath,testbase+'*')))
testlen  = len(testdirs)

# allocate empty arrays and dictionary for convergence results
results = np.zeros([testlen, 12])
summary = {}
summary['frame'] = frame

# =====================================================================================

print '\n=============================================\n'
print 'looking at frame ', str(frame)

create_dir(savepath)
create_dir(figspath)
create_dir(repopath)

# Convergence test at interest_frame
for m,enddir in enumerate(xrange(7,13,1)):
    dirs = os.path.join(testpath,testbase+str(enddir))
    print '\n-------------------------------------------'
    print '\ndir ', dirs+''

    # load first pyclaw solution q
    solution = Solution()
    solution.read(frame,path=dirs,file_format='petsc',read_aux=False)

    # self convergence test
    if m<testlen-1:
        solution_test = Solution()
        dirs2 = os.path.join(testpath,testbase+str(enddir+1))
        print dirs2
        solution_test.read(frame,path=dirs2,file_format='petsc',read_aux=False)
        self_difference = self_convergence(solution,solution_test,debug=debug_flag)

    # write results to compiling array of results (convinience array to save in matlab)
    results[m,0]     = solution.patch.num_cells_global[0]
    results[m,1]     = np.prod(solution.state.grid.delta)
    results[m,2:4]   = self_difference

# pth order
results[0:-1,4] = np.log2(results[0:-1,3]/results[1:,3])
headers = ['$n_{cells}$','$h$','$E_s(h)$','$p_s$']
results = results[:,[0,1,3,4]]
try:
    from tabulate import tabulate
    strt = tabulate(results,headers=headers,tablefmt="latex",floatfmt="4.3e")
    dstfile = os.path.join(repopath,'table_errors.tex')
    if debug:
        print strt
except:
    pass

rate = np.zeros([2,2])
x = results[np.where(results[:,3]>=1.5),0]
r = results[np.where(results[:,3]>=1.5),2].T
rate[0,:] = linear_fit(np.log10(x),np.log10(r))

x = results[np.where(results[:,3]>=1.5),1]
r = results[np.where(results[:,3]>=1.5),2].T
rate[1,:] = linear_fit(np.log10(x),np.log10(r))

strt = strt + '\n'+ tabulate(rate,headers=['$p$','$c$'],tablefmt="latex",floatfmt="4.3e")
f = open(dstfile,'w`')
f.write(strt)
f.close()

summary['convergence'] = results
summary['log_convergence'] = np.log10(results)
summary['linear_fit'] = rate

np.save(os.path.join(savepath,outprefix+'convergence_f'+str(frame)),results)
savemat(os.path.join(savepath,outprefix+'summary_f'+str(frame)),summary)

# plot results
convergence_cases = ['self']
save_res_name = 'convergence_'

x = results[:,0]
for n,case in enumerate(convergence_cases):
    y = results[:,2*(n+1)]
    plot_convergence(x,y,rate=rate[n,:],log=True,outdir=figspath,figname=save_res_name,case=case,tolatex=True)
