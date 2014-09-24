import os
from glob import glob

from clawpack.pyclaw import Solution

from scipy.io import loadmat,savemat
import numpy as np
import matplotlib
# set matplotlib to work over X-forwarding
matplotlib.use('Agg')
from matplotlib import pylab as plt

def debug_convergence_message(results,restype):
    print '\n'+restype+' (absolute)\t\t'+str(results[0])
    print restype+' (relative)\t\t'+str(results[1])

def self_convergence(claw_old,claw_new,q=0,debug=False):
    pfinal1 = claw_old.state.get_q_global()[q]
    pfinal2 = claw_new.state.get_q_global()[q]

    delta = claw_new.state.grid.delta

    # average claw new
    if claw_new.state.num_dim==1:
        pfinal2 = (pfinal2[::2] + pfinal2[1::2])/2.0
    elif claw_new.state.num_dim==2:
        pfinal2 = (pfinal2[::2,::2] + pfinal2[1::2,::2] + pfinal2[::2,:1:2] + pfinal2[1::2,1::2])/4.0
    elif claw_new.state.num_dim==3:
        pfinal2 = (pfinal2[::2,::2,::2] + pfinal2[1::2,::2,::2] + pfinal2[::2,:1:2,::2] + pfinal2[1::2,1::2,::2] + 
            pfinal2[::2,::2,1::2] + pfinal2[1::2,::2,1::2] + pfinal2[::2,:1:2,1::2] + pfinal2[1::2,1::2,1::2])/8.0

    pfinal2 = pfinal2.reshape(-1)
    pfinal1 = pfinal1.reshape(-1)

    self_difference_absolute = np.prod(delta)*np.linalg.norm(pfinal1-pfinal2,ord=1)
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

def plot_analytic_convergence(x1,q1,x2,q2,outname='outplot',outdir='./figures',label1='exact',label2='numeric',format='eps',titleon=True):

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

def plot_convergence(x,y,log=True,outdir='./',figname='_plot',case='',title='',format='eps',tolatex=False):
    create_dir(outdir)
    log_name  = ''
    log_label = ''
    log_title = ''
    if not tolatex:
        log_title = 'Log Log '
    if log:
        log_name  = 'loglog_'
    
    plt.figure()
    if log:
        plt.loglog(x,y,'o:')
    else:
        plt.plot(x,y,'o:')
    if not tolatex:
        plt.title(title+log_title+case+' convergence')
    
    plt.xlabel('$'+log_label+'n_{cells}$')
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
    dstfile = os.path.join(repopath,'summary_'+restype+'_figures.tex')
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

def print_summary(rate,cases,restype='relative',tolatex=False):

    strt = ''
    for m,case in enumerate(cases):
        strt = strt + case +'\t' + str(np.abs(rate[m,0]))

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
        for m,case in enumerate(cases):
            strt = strt + '\n' + '\t'+case+'\t&\t' + str(np.abs(rate[m,0])) + r' \\'
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
    A = np.vstack([x,np.ones(len(x))]).T
    m,c = np.linalg.lstsq(A,y)[0]
    return [m,c]

def zoom_max(q,x,delta):
    indmax = np.where(q==q.max())[0][0]
    zoom_range = [x[indmax]-delta, x[indmax]+delta]
    return indmax,zoom_range

def convergence(suffix,savepath='./',figpath='./',frame=0,outdir='./',basefile='_output_',selfc=True,analytic=False,point=False,genmat=False,matfile=None,delta=2.0,debug=True,q=1):

    summary = {}

    results = np.zeros([len(suffix),2+2*np.sum([selfc,analytic,point])])
    
    outs = ''
    
    if selfc: outs = outs+'.self'
    
    if analytic:
        outs = outs + '.analytic'
        print '\nloading mat file ', matsrc
        matdict = loadmat(os.path.join(matpath,matsrc),variable_names=['Q','X','dx'])

        q_exact = matdict['Q'][:,-1]
        x_exact = matdict['X'][:,-1]

        indmax,zoom_range = zoom_max(q_exact,x_exact,delta)

        if debug:
            print indmax, q_exact[indmax], q_exact.max()
            print indmax, x_exact[indmax]
            print zoom_range    
    
    if point: outs = outs+'.point'

    for m,s in enumerate(suffix):
        srcdir = os.path.join(outdir,basefile+str(s))
        print '\n-------------------------------------------'
        print '\ndir ', srcdir+''

        # load first pyclaw solution q
        solution = Solution()
        solution.read(frame,path=srcdir,file_format='petsc',read_aux=False)

        results[m,0] = np.prod(solution.patch.num_cells_global)
        results[m,1] = np.prod(solution.state.grid.delta)

        if selfc:
            if m<(len(suffix)-1):
                solution_test = Solution()
                srcdir = os.path.join(outdir,basefile+str((s+1)))
                solution_test.read(frame,path=srcdir,file_format='petsc',read_aux=False)
                self_difference = self_convergence(solution,solution_test,debug=debug,q=q)
                results[m,2:4] = self_difference

        if analytic:
            global_difference = global_analytic_convergence(solution,q_exact,x_exact,q=q,debug=debug,frame=frame,
                outdir=figpath,plot=True)
            results[m,4:6] = global_difference

        if point:
            point_difference  = point_analytic_convergence(solution, q_exact,x_exact,q=q,debug=debug)
            results[m,6:8] = point_difference

        if genmat:
            summary['frame'] = frame
            summary['t'+str(m)] = solution.state.t
            summary['Q'+str(m) ] = solution.state.get_q_global()[q]
            summary['X'+str(m) ] = solution.state.grid.x.centers
            summary['dX'+str(m)] = solution.state.grid.x.delta
            summary['N'+str(m) ] = solution.state.grid.x.num_cells

    outs = outs.split('.')[1::]

    return results,summary,outs

def set_env(savepath):
    create_dir(savepath)
    figs = os.path.join(savepath,'_figures')
    reps = os.path.join(savepath,'_reports')
    create_dir(figs)
    create_dir(reps)

    return figs,reps
# def summarize(results,figpath,reppath,rate_calc_range=None,summary={},genmat=False):

# =====================================================================================
# matlab source path and filename
matpath = '/simdesk/sandbox/emclaw/analysis/results/2D'
matsrc  = 'analytic_centers_all_hd_exact_32768.mat'

# pyclaw source path and filenames
testpath = '/media/shaheen/dev/results/maxwell_2d'
testbase = '_output_'

# output peraphernalia
savepath = '/simdesk/sandbox/emclaw/analysis/results/2D/test'
outprefix = 'f55_matgrid_'

# define frame of interest for analytic vs numeric convergence
frame = 55
qn = 1

frange = [7,13]
suffix = range(frange[0],frange[1]+1,1)

genmat = False

figpath,reppath = set_env(savepath)

results,summary,cases = convergence(suffix,frame=frame,outdir=testpath,basefile=testbase,q=qn,savepath=savepath,figpath=figpath)
repopath = reppath
x = results[:,0].transpose()
r = results[:,3::2].transpose()

x = x[2:(len(suffix)-1)]
r = r[:,2:(len(suffix)-1)]

rate = np.zeros([len(r),2])
print x
print r
print r.shape
for k in xrange(0,len(r),1):
    print k
    rate[k,:] = linear_fit(np.log10(x),np.log10(r[k]))

summary['convergence'] = results
summary['log_convergence'] = np.log10(results)
summary['linear_fit'] = rate

np.save(os.path.join(savepath,outprefix+'convergence_f'+str(frame)),results)
savemat(os.path.join(savepath,outprefix+'summary_f'+str(frame)),summary)

print rate,cases

#print_summary(rate,cases,restype='relative',tolatex=True)

restype = 'relative'
save_res_name = 'convergence'+'_'+restype

x = results[:,0]

if restype=='absolute':
    results = results[:,2::2]
else:
    results = results[:,3::2]

for n,case in enumerate(cases):
    y = results[:,n]
    plot_convergence(x,y,log=True,outdir=figpath,figname=save_res_name,case=case,tolatex=True)
