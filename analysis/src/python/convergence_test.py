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

# matlab source path and filename
matpath = '/simdesk/sandbox/emclaw/results/1D/convergence/summary'
matsrc  = 'analytic_centers_all_hd_exact_131072.mat'

# pyclaw source path and filenames
testpath = '/simdesk/sandbox/emclaw/results/1D/convergence_averaged/'
testbase = '_output_'

# output peraphernalia
savepath = '/simdesk/sandbox/emclaw/results/1D/convergence_averaged/summary'
figspath = '/simdesk/sandbox/emclaw/results/1D/convergence_averaged/summary/figures'
repopath = '/simdesk/sandbox/emclaw/results/1D/convergence_averaged/summary/reports'
outprefix = 'f55_matgrid_'

# define frame of interest for analytic vs numeric convergence
frame = 55
qn = 0

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

    pfinal2 = pfinal2[0]
    pfinal1 = pfinal1[0]

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

def local_analytic_convergence(claw,qmat,xmat,plot=False,local_range=[],frame=50,q=0,outdir='./figures',debug=False):
    qclaw  = claw.state.get_q_global()[q]
    xclaw  = claw.state.grid.x.centers

    qclaw = reduce_dim(qclaw)

    p = qmat.shape[0]/qclaw.shape[0]

    qmat = qmat[::p]
    xmat = xmat[::p]

    # take a portion of the array defined by the limits in local_range
    qmat  = qmat[(xclaw>=local_range[0])*(xclaw<=local_range[1])]
    qclaw = qclaw[(xclaw>=local_range[0])*(xclaw<=local_range[1])]
    xmat = xmat[(xclaw>=local_range[0])*(xclaw<=local_range[1])]
    xclaw  = xclaw[(xclaw>=local_range[0])*(xclaw<=local_range[1])]

    # reshape, this might be unnecesary in 1D
    qmat  = qmat.reshape(-1)
    xmat  = xmat.reshape(-1)
    qclaw = qclaw.reshape(-1)
    xclaw = xclaw.reshape(-1)

    # calculate the difference
    local_analytic_difference_absolute = np.prod(claw.state.grid.x.delta)*np.linalg.norm(qclaw-qmat,ord=1)
    local_analytic_difference_relative = local_analytic_difference_absolute/np.linalg.norm(qmat,ord=1)

    local_analytic_difference = [local_analytic_difference_absolute, local_analytic_difference_relative
    ]
    if plot:
        figname = 'analytics_numeric_local_f'+str(frame)+'nc_'+str(int(claw.state.grid.x.num_cells))
        plot_analytic_convergence(xmat,qmat,xclaw,qclaw,outname=figname,outdir=outdir,titleon=False)

    if debug:
        debug_convergence_message(local_analytic_difference,'analytic local')

    return local_analytic_difference

def global_analytic_convergence(claw,qmat,xmat,plot=False,frame=50,q=0,outdir='./figures',debug=False):
    qclaw  = claw.state.get_q_global()[q]
    xclaw  = claw.state.grid.x.centers

    qclaw = reduce_dim(qclaw)

    p = qmat.shape[0]/qclaw.shape[0]
    
    qmat = qmat[::p]
    xmat = xmat[::p]
    
    qclaw = qclaw.reshape(-1)
    xlcaw = xclaw.reshape(-1)
    qmat = qmat.reshape(-1)
    xmat = xmat.reshape(-1)

    global_analytic_difference_absolute = np.prod(claw.state.grid.x.delta)*np.linalg.norm(qclaw-qmat,ord=1)
    global_analytic_difference_relative = global_analytic_difference_absolute/np.linalg.norm(qmat,ord=1)

    global_analytic_difference = [global_analytic_difference_absolute,global_analytic_difference_relative]

    if plot:
        figname = 'analytics_numeric_global_f'+str(frame)+'nc_'+str(int(claw.state.grid.x.num_cells))
        plot_analytic_convergence(xmat,qmat,xclaw,qclaw,outname=figname,outdir=outdir)

    if debug:
        debug_convergence_message(global_analytic_difference,'analytic global')
        
    return global_analytic_difference

def point_analytic_convergence(claw,qmat,xmat,q=0,debug=False):
    qclaw  = claw.state.get_q_global()[q]
    xclaw  = claw.state.grid.x.centers

    qclaw = reduce_dim(qclaw)

    p = qmat.shape[0]/qclaw.shape[0]
    
    qmat = qmat[::p]
    xmat = xmat[::p]

    index = np.where(qmat>=qmat.max())[0]

    point_analytic_difference_absolute = np.abs(qclaw[index] - qmat[index])[0]
    point_analytic_difference_relative = point_analytic_difference_absolute/qmat[index][0]
    
    point_analytic_difference = [point_analytic_difference_absolute,point_analytic_difference_relative]
    print point_analytic_difference
    if debug:
        debug_convergence_message(point_analytic_difference,'analytic point')


        print '\nPoint wise comparisson (static)'
        print '\n                Claw         Analytic'
        print '------------------------------------------'
        print 'index max   ', index[0], '        ',index[0]
        print 'q compare   ', qclaw[index][0], q_exact[index][0]
        print 'x compare   ', xclaw[index][0], x_exact[index][0]
        print '------------------------------------------'
        print 'x diff      ', np.abs(xclaw[index] - x_exact[index])[0]
        print 'q diff      ', np.abs(qclaw[index] - q_exact[index])[0],'\n'
    
    return point_analytic_difference

def point_seek_analytic_convergence(claw,qmat,xmat,seek_index,q=0,debug=False):
    qclaw  = claw.state.get_q_global()[q]
    xclaw  = claw.state.grid.x.centers

    qclaw = reduce_dim(qclaw)
    
    indclaw = np.where(xclaw>=x_exact[seek_index])[0]

    if np.shape(indclaw)[0]>1:
        indclaw = indclaw[0]

    point_seek_difference_absolute = np.abs(qclaw[indclaw] - q_exact[indmax])
    point_seek_difference_relative = point_seek_difference_absolute/q_exact[indmax]

    point_seek_difference = [point_seek_difference_absolute,point_seek_difference_relative]

    if debug:
        debug_convergence_message(point_seek_difference,'analytic seek')

        print '\nPoint wise comparisson (seek)'
        print '\n                Claw         Analytic'
        print '------------------------------------------'
        print 'index max   ', indclaw, '        ',seek_index
        print 'q compare   ', qclaw[indclaw], q_exact[seek_index]
        print 'x compare   ', xclaw[indclaw], x_exact[seek_index]
        print '------------------------------------------'
        print 'x diff      ', np.abs(xclaw[indclaw] - x_exact[seek_index])
        print 'q diff      ', np.abs(qclaw[indclaw] - q_exact[seek_index]),'\n'

    return point_seek_difference

def plot_analytic_convergence(x1,q1,x2,q2,outname='outplot',outdir='./figures',label1='exact',label2='numeric',format='eps',titleon=False):

    create_dir(outdir)

    plt.figure()
    plt.plot(x1,q1,':',label=label1)
    plt.plot(x2,q2,':',label=label2)

    if titleon: plt.title('Analytic vs Numeric comparisson at $n_c = '+str(q1.shape[0])+'$')
    plt.xlabel('$x$ (a.u.)')
    plt.ylabel('$Q$ (a.u.)')
    plt.legend(frameon=False)

    plt.draw()
    plt.savefig(os.path.join(outdir,outname+'.'+format),format=format,dpi=320,bbox_inches='tight')
    plt.close()

    return

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

def print_summary(rate,restype='relative',tolatex=False):

    strt = ''
    strt = strt + '\n' + 'self convergence  ' + str(np.abs(rate[3,0]))
    strt = strt + '\n' + 'global analytic   ' + str(np.abs(rate[2,0]))
    strt = strt + '\n' + 'local analytic    ' + str(np.abs(rate[1,0]))
    strt = strt + '\n' + 'point (static)    ' + str(np.abs(rate[0,0]))
    strt = strt + '\n' + 'point (seek)      ' + str(np.abs(rate[4,0]))

    print '\n============================================='
    print 'Summary ('+restype+')'
    print '\n============================================='
    print  strt
    print '\============================================='

    if tolatex:
        strt = r'\begin{table}'
        strt = strt + '\n' + r'\centering'
        strt = strt + '\n' + r'\begin{tabular}{ll}'
        strt = strt + '\n' + '\t\\hline'
        strt = strt + '\n' + '\t\\emph{Convergence measurement} \t&\t \\emph{Rate}' + r' \\'
        strt = strt + '\n' + '\t\\hline'
        strt = strt + '\n' + '\tself convergence  \t&\t' + str(np.abs(rate[3,0])) + r' \\'
        strt = strt + '\n' + '\tglobal analytic   \t&\t' + str(np.abs(rate[2,0])) + r' \\'
        strt = strt + '\n' + '\tlocal analytic    \t&\t' + str(np.abs(rate[1,0])) + r' \\'
        strt = strt + '\n' + '\tpoint (static)    \t&\t' + str(np.abs(rate[0,0])) + r' \\'
        strt = strt + '\n' + '\tpoint (seek)      \t&\t' + str(np.abs(rate[4,0])) + r' \\'
        strt = strt + '\n' + '\t\\hline'
        strt = strt + '\n' + r'\end{tabular}'
        strt = strt + '\n' + r'\caption{Summary, convergence results ($\Delta q/q_{exact}$)}'
        strt = strt + '\n' +r'\end{table}'
        create_dir(repopath)
        dstfile = os.path.join(repopath,'summary_'+restype+'_table.tex')
        f = open(dstfile,'w')
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

# load matlab file for analytic convergence
print '\nloading mat file ', matsrc
matdict = loadmat(os.path.join(matpath,matsrc),variable_names=['Q','X','dx'])

q_exact = matdict['Q'][:,-1]
x_exact = matdict['X'][:,-1]

# find location of Q_exact max and set local_range
indmax = np.where(q_exact==q_exact.max())[0][0]
zoom_range = [x_exact[indmax]-local_range_delta, x_exact[indmax]+local_range_delta]
if debug_flag:
    print indmax, q_exact[indmax], q_exact.max()
    print indmax, x_exact[indmax]
    print zoom_range

# =====================================================================================

print '\n=============================================\n'
print 'looking at frame ', str(frame)

create_dir(savepath)
create_dir(figspath)
create_dir(repopath)

# Convergence test at interest_frame
for m,enddir in enumerate(xrange(7,15,1)):
    dirs = os.path.join(testpath,testbase+str(enddir))
    print '\n-------------------------------------------'
    print '\ndir ', dirs+''

    # load first pyclaw solution q
    solution = Solution()
    solution.read(frame,path=dirs,file_format='petsc',read_aux=False)

    # get the point, local and global analytic convergence
    local_difference  = local_analytic_convergence(solution, q_exact,x_exact,q=qn,debug=debug_flag,frame=frame,outdir=figspath,local_range=zoom_range,plot=True)
    global_difference = global_analytic_convergence(solution,q_exact,x_exact,q=qn,debug=debug_flag,frame=frame,outdir=figspath,plot=True)

    # self convergence test
    if m<testlen-1:
        solution_test = Solution()
        dirs2 = os.path.join(testpath,testbase+str(enddir+1))
        print dirs2
        solution_test.read(frame,path=dirs2,file_format='petsc',read_aux=False)
        self_difference = self_convergence(solution,solution_test,debug=debug_flag,dim='x')

    # write results to compiling array of results (convinience array to save in matlab)
    results[m,0]     = solution.patch.num_cells_global[0]
    results[m,1]     = solution.state.grid.delta[0]
    results[m,2:4]   = local_difference
    results[m,4:6]   = global_difference
    results[m,6:8]   = self_difference

    summary['t'+str(m)] = solution.state.t

    # build a dictionary with Q at frame, for postprocessing
    summary['Q'+str(m) ] = solution.state.get_q_global()
    summary['aux'+str(m)] = solution.state.get_aux_global()
    summary['X'+str(m) ] = solution.state.grid.x.centers
    summary['dX'+str(m)] = solution.state.grid.x.delta
    summary['N'+str(m) ] = solution.state.grid.x.num_cells

# pth order
results[0:-1,8:11] = np.log2(results[0:-1,[3,5,7]]/results[1:,[3,5,7]])
headers = ['$n_{cells}$','$h$','$E_e(h)$','$p_e$','$E_s(h)$','$p_s$']
results = results[:,[0,1,5,9,7,10]]
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

x = results[np.where(results[:,5]>=1.5),0]
r = results[np.where(results[:,5]>=1.5),4].T
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
convergence_cases = ['global','self']
save_res_name = 'convergence_'

x = results[:,0]
for n,case in enumerate(convergence_cases):
    y = results[:,2*(n+1)]
    plot_convergence(x,y,rate=rate[n,:],log=True,outdir=figspath,figname=save_res_name,case=case,tolatex=True)
