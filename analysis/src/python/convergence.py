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

class Convergence(object):
        
    def debug_message(self,results,restype):
        print '\n'+restype+' (absolute)\t\t'+str(results[0])
        print restype+' (relative)\t\t'+str(results[1])

    def reduce_dim(self,q):
        if len(q.shape)==2:
            q = q[:,q.shape[1]/2]
        elif len(q.shape)==3:
            q = q[:,q.shape[1]/2,q.shape[2]/2]

        return q

    def average_reduce(self,q,p=0):

        for i in range(int(p)):
            if self.num_dim==1:
                q = (q[::2] + q[1::2])/2.0
            
            if self.num_dim==2:
                q = (q[::2,::2] + q[1::2,::2] + q[::2,:1:2] + q[1::2,1::2])/4.0
            
            if self.num_dim==3:
                q = (q[::2,::2,::2] + q[1::2,::2,::2] + q[::2,:1:2,::2] + q[1::2,1::2,::2] + 
                    q[::2,::2,1::2] + q[1::2,::2,1::2] + q[::2,:1:2,1::2] + q[1::2,1::2,1::2])/8.0
        return q

    def errors(self,qold,qnew,delta):
        difference = np.zeros([2])
        difference[0] = np.prod(delta)*np.linalg.norm((qold - qnew),ord=1)
        difference[1] = difference[0]/np.linalg.norm(qnew,ord=1)
        return difference

    def self_convergence(self,qold,qnew,delta):

        p = np.log2(qnew.shape[0]/qold.shape[0])
        qnew = self.average_reduce(qnew,p)

        qnew = qnew.reshape(-1)
        qold = qold.reshape(-1)

        self_difference = self.errors(qnew,qold,delta)

        if self.debug:
            self.debug_message(self_difference,'self')
        
        return self_difference

    def analytic_convergence(self,xclaw,qclaw,xmat,qmat,delta,local=False,finest=False):
        affix = 'global'
        local_difference = np.zeros([2])
        p = np.log2(qmat.shape[0]/qclaw.shape[0])

        qmat = self.average_reduce(qmat,p)
        xmat = self.average_reduce(xmat,p)
        # take a portion of the array defined by the limits in local_range
        if local:
            qmat  = qmat[(xclaw>=self.zoom_range[0])*(xclaw<=self.zoom_range[1])]
            qclaw = qclaw[(xclaw>=self.zoom_range[0])*(xclaw<=self.zoom_range[1])]
            xmat  = xmat[(xclaw>=self.zoom_range[0])*(xclaw<=self.zoom_range[1])]
            xclaw = xclaw[(xclaw>=self.zoom_range[0])*(xclaw<=self.zoom_range[1])]
            affix = 'local'
        if finest:
            affix = 'finest_'+affix

        # reshape, this might be unnecesary in 1D
        qmat  = qmat.reshape(-1)
        xmat  = xmat.reshape(-1)
        qclaw = qclaw.reshape(-1)
        xclaw = xclaw.reshape(-1)

        # calculate the difference
        analytic_difference = self.errors(qclaw,qmat,delta)

        if self.plot:
            figname = '_'+affix+'_analytic_'+str(self.frame)+'_p_'+str(int(p))
            self.plot_solutions(xmat,qmat,xclaw,qclaw,figname=figname,titleon=False)

        if self.debug:
            self.debug_message(analytic_difference,'analytic '+affix)

        return analytic_difference

    def plot_solutions(self,x1,q1,x2,q2,figname='outplot',label1='exact',label2='numeric',titleon=False):

        self.create_dir(self.plotdir)

        plt.figure()
        plt.plot(x1,q1,':',label=label1)
        plt.plot(x2,q2,':',label=label2)

        if titleon: plt.title('Analytic vs Numeric comparisson at $n_c = '+str(q1.shape[0])+'$')
        plt.xlabel('$x$ (a.u.)')
        plt.ylabel('$Q$ (a.u.)')
        plt.legend(frameon=False)

        plt.draw()
        plt.savefig(os.path.join(self.plotdir,figname+'.'+self.plot_format),format=self.plot_format,dpi=320,bbox_inches='tight')
        plt.close()

        return

    def create_dir(self,directory):
        if not os.path.exists(directory): os.makedirs(directory)
        return

    def plot_convergence(self,x,y,rate=None,log=True,figname='_plot',case='',title='',tolatex=False):
        self.create_dir(self.plotdir)
        log_name  = ''
        log_label = ''
        log_title = ''
        if not tolatex:
            log_title = 'Log Log '
        if log:
            log_name  = 'loglog_'

        plt.figure()
        if rate is not None:
            p = self.p_line_range
            m = rate[0]
            c = rate[1]
            plt.hold(True)
            plt.loglog(x[p[0]:p[1]],10.0**(m*np.log10(x[p[0]:p[1]])+c-2.0),'r')
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
        fig_save_name = figname+'_'+log_name+case+'.'+self.plot_format
        figpath = os.path.join(self.plotdir,fig_save_name)
        plt.savefig(figpath,format=self.plot_format,dpi=320,bbox_inches='tight')
        plt.close()

        if tolatex:
            caption = ''
            if log:
                caption = 'Log Log '
            caption = caption + 'plot of '+case.replace("_"," ")+' convergence test'
            caption.capitalize()
            self.gen_latex_fig(figpath,caption=caption)
        return

    def gen_latex_fig(self,figpath,caption=''):
        self.create_dir(self.reportdir)
        dstfile = os.path.join(self.reportdir,'summary_figures.tex')
        f = open(dstfile,'a')
        strt = r'\begin{figure}[h!]'
        strt = strt + '\n \t' + r'\centering'
        strt = strt + '\n \t' + r'\includegraphics[width=0.75\textwidth]{'+os.path.relpath(figpath,self.reportdir)+'}'
        strt = strt + '\n \t' + r'\caption{'+caption+'}'
        strt = strt + '\n \t' + r'\label{fig:xxxx}'
        strt = strt + '\n' + r'\end{figure}'
        strt = strt + '\n\n'
        f.write(strt)
        f.close()
        return

    def linear_fit(self,x,y):
        A = np.vstack([x,np.ones(x.size)]).T
        m,c = np.linalg.lstsq(A,y)[0]
        return m,c

    def getdirs(self):
        testdirs = sorted(glob(os.path.join(self.testdir,self.basedir+'*')))

        if self.debug:
            print 'Test dirs are:'
            print testdirs

        return testdirs

    def update_dir(self):
        setattr(self, 'plotdir', os.path.join(self.savedir,'_plots'))
        setattr(self, 'reportdir', os.path.join(self.savedir,'_report'))

    def getClaw(self,dirs,frame=None,qn=None,homogeneous=None):
        if frame is None:
            frame = self.frame
        if qn is None:
            qn = self.qn
        if homogeneous is None:
            homogeneous = self.homogeneous

        solution = Solution()

        solution.read(frame,path=dirs,file_format=self.file_format,read_aux=self.homogeneous)

        qclaw  = solution.state.get_q_global()[qn]
        xclaw  = solution.state.grid.x.centers
        delta  = solution.state.grid.delta

        if homogeneous:
            print 'dividing by aux'
            qclaw = qclaw/(solution.state.get_aux_global()[qn])

        return qclaw,xclaw,delta

    def getAux(self,dirs,frame=None,qn=None):
        if frame is None:
            frame = self.frame
        if qn is None:
            qn = self.qn

        solution = Solution()
        solution.read(frame,path=dirs,file_format=self.file_format,read_aux=True)

        auxclaw  = solution.state.get_aux_global()[qn]
        xclaw    = solution.state.grid.x.centers
        delta    = solution.state.grid.delta

        return auxclaw,xclaw,delta

    def load_matlab(self):

        matdict = loadmat(self.matsrc,variable_names=self.matvar)

        for var in self.matvar:
            setattr(self, var, matdict[var][:,-1])

        q_exact = matdict['qs'][:,-1]
        x_exact = matdict['xs'][:,-1]

        return q_exact,x_exact

    def compare(self,dir1,dir2,frame=None,plot=False,affix='',aux=False):
        if frame is None:
            frame = self.frame

        if aux:
            q1,x1,d1 = self.getAux(dir1,frame)
            q2,x2,d2 = self.getAux(dir2,frame)
        else:
            q1,x1,d1 = self.getClaw(dir1,frame)
            q2,x2,d2 = self.getClaw(dir2,frame)

        compare = self.self_convergence(q1,q2,d1)

        if plot:
            self.plot_solutions(x1,q1,x2,q2,figname=affix+'frame_'+str(frame),label1='s1',label2='s2')

        return compare

    def makedirs(self):
        self.create_dir(self.savedir)
        self.create_dir(self.plotdir)
        self.create_dir(self.reportdir)

    def __init__(self,testdir='./',basedir='_output_',savedir=None,frame=0):
        self.testdir = testdir
        self.basedir = basedir
        self.compare_dir = testdir
        self.compare = False
        self.frame   = frame
        self.basemin = 7
        self.basemax = 15
        self.summary = {}
        self.debug   = False
        self.plot    = True
        if savedir is None:
            savedir    = testdir
        self.savedir   = savedir
        self.plotdir   = None
        self.reportdir = None
        self.update_dir()
        self.qn           = 0
        self.pth          = 1.5
        self.plot_format  = 'eps'
        self.p_line_range = [2,7]
        self.file_format  = 'petsc'
        self.homogeneous  = False

class Errors1D(Convergence):
    def __init__(self,testdir='./',basedir='_output_',savedir=None,frame=0):
        self.num_dim = 1
        self.matsrc  = './matlab.mat'
        self.matvar  = ['qs','xs']
        self.finesrc = './_output_finest'
        self.range_delta = 2.0
        self.convergence_cases = ['global','finest','self']
        super(Errors1D,self).__init__(testdir,basedir,savedir,frame)

    def local_range(self,q_exact,x_exact):
        indmax = np.where(q_exact==q_exact.max())[0][0]
        zoom_range = [x_exact[indmax]-self.range_delta, x_exact[indmax]+self.range_delta]

        setattr(self, 'zoom_range', zoom_range)
        setattr(self, 'indmax', indmax)

        if self.debug:
            print indmax, q_exact[indmax], q_exact.max()
            print indmax, x_exact[indmax]
            print zoom_range

        return indmax,zoom_range

    def convergence(self):
        self.update_dir()

        testlen = self.basemax-self.basemin+1
        if self.compare:
            add = 1
        else:
            add = 0

        results = np.zeros([testlen,12+add])

        qmat,xmat   = self.load_matlab()
        indmax,zoom_range = self.local_range(qmat,xmat)

        self.makedirs()

        # load solution at finest grid to be used as reference
        qfinest,xfinest,deltafinest = self.getClaw(self.finesrc)

        for m,enddir in enumerate(xrange(self.basemin,self.basemax+1,1)):
            dirs = os.path.join(self.testdir,self.basedir+str(enddir))

            print '\n-------------------------------------------'
            print '\ndir ', dirs+''

            # load first pyclaw solution qn
            qclaw,xclaw,delta = self.getClaw(dirs)

            # get the error with respect to the exact solution (local and global)
            local_difference  = self.analytic_convergence(xclaw,qclaw,xmat,qmat,delta,local=True)
            global_difference = self.analytic_convergence(xclaw,qclaw,xmat,qmat,delta,local=False)

            # get the error with respect to the finest solution
            finest_difference = self.self_convergence(qclaw,qfinest,delta)
            
            # compare against another claw solution
            if self.compare:
                dirs2 = os.path.join(self.compare_dir,self.basedir+str(enddir))
                qclaw2,xclaw2,delta2 = self.getClaw(dirs2,homogeneous=False)
                compare_difference = self.errors(qclaw,qclaw2,delta)
                results[m,-1] = compare_difference[1]
                figname = '_compare_'+str(self.frame)+'_m_'+str(int(m))
                self.plot_solutions(xclaw,qclaw,xclaw2,qclaw2,figname=figname,titleon=False)

            # get the error with respect to the refined solution
            if enddir<self.basemax:
                dir_refined = os.path.join(self.testdir,self.basedir+str(enddir+1))
                print '\n-------------------------------------------'
                print '\ndir ', dir_refined

                qref,xref,deltaref = self.getClaw(dir_refined)
                sfiner_difference  = self.self_convergence(qclaw,qref,delta)

            # write results to compiling-array
            results[m,0]     = 2**(enddir)
            results[m,1]     = np.prod(delta)
            results[m,2:4]   = global_difference
            results[m,4:6]   = finest_difference
            results[m,6:8]   = sfiner_difference

        # calc pth order
        results[0:-1,8:11] = np.log2(results[0:-1,[3,5,7]]/results[1:,[3,5,7]])
        if self.compare:
            headers = ['$n_{cells}$','$h$','$E_e(h)$','$p_e$','$E_f(h)$','$p_f$','$E_s(h)$','$p_s$','comp']
            results = results[:,[0,1,3,8,5,9,7,10,12]]
        else:
            headers = ['$n_{cells}$','$h$','$E_e(h)$','$p_e$','$E_f(h)$','$p_f$','$E_s(h)$','$p_s$']
            results = results[:,[0,1,3,8,5,9,7,10]]

        try:
            from tabulate import tabulate
            strt = tabulate(results,headers=headers,tablefmt="latex",floatfmt="4.3e")
            dstfile = os.path.join(self.reportdir,'table_errors.tex')
            if self.debug:
                print strt
        except:
            print 'Unable to create table of errors, have you done pip install tabulate'
            pass

        rate = np.zeros([3,2])
        x = results[np.where(results[:,3]>=self.pth),0]
        r = results[np.where(results[:,3]>=self.pth),2].T
        rate[0,:] = self.linear_fit(np.log10(x),np.log10(r))

        x = results[np.where(results[:,5]>=self.pth),0]
        r = results[np.where(results[:,5]>=self.pth),4].T
        rate[1,:] = self.linear_fit(np.log10(x),np.log10(r))

        x = results[np.where(results[:,7]>=self.pth),0]
        r = results[np.where(results[:,7]>=self.pth),6].T
        rate[2,:] = self.linear_fit(np.log10(x),np.log10(r))

        try:
            strt = strt + '\n'+ tabulate(rate,headers=['$p$','$c$'],tablefmt="latex",floatfmt="4.3e")
            f = open(dstfile,'w`')
            f.write(strt)
            f.close()
            if self.debug:
                print strt
        except:
            pass

        self.summary['convergence'] = results
        self.summary['log_convergence'] = np.log10(results)
        self.summary['linear_fit'] = rate

        np.save(os.path.join(self.savedir,'convergence_f'+str(self.frame)),results)
        savemat(os.path.join(self.savedir,'summary_f'+str(self.frame)),self.summary)

        # plot results
        if self.plot:
            save_res_name = '_convergence_'
            x = results[:,0]
            for n,case in enumerate(self.convergence_cases):
                y = results[:,2*(n+1)]
                self.plot_convergence(x,y,rate=rate[n,:],log=True,figname=save_res_name,case=case,tolatex=True)

class Errors2D(Convergence):
    def __init__(self,testdir='./',basedir='_output_',savedir=None,frame=0):
        self.num_dim = 2
        self.finesrc = './_output_finest'
        self.convergence_cases = ['finest','self']
        super(Errors2D,self).__init__(testdir,basedir,savedir,frame)
        self.qn = 1

    def convergence(self):
        self.update_dir()

        testlen = self.basemax-self.basemin+1

        results  = np.zeros([testlen,12])

        self.makedirs()

        # load solution at finest grid to be used as reference
        qfinest,xfinest,deltafinest = self.getClaw(self.finesrc)

        for m,enddir in enumerate(xrange(self.basemin,self.basemax+1,1)):
            dirs = os.path.join(self.testdir,self.basedir+str(enddir))

            print '\n-------------------------------------------'
            print '\ndir ', dirs+''

            # load first pyclaw solution qn
            qclaw,xclaw,delta = self.getClaw(dirs)

            # get the error with respect to the finest solution
            finest_difference = self.self_convergence(qclaw,qfinest,delta)
            
            # get the error with respect to the refined solution
            if enddir<self.basemax:
                dir_refined = os.path.join(self.testdir,self.basedir+str(enddir+1))
                print '\n-------------------------------------------'
                print '\ndir ', dir_refined

                qref,xref,deltaref = self.getClaw(dir_refined)
                sfiner_difference  = self.self_convergence(qclaw,qref,delta)

            # write results to compiling-array
            results[m,0]     = 2**(enddir)
            results[m,1]     = np.prod(delta)
            results[m,2:4]   = finest_difference
            results[m,4:6]   = sfiner_difference

        # calc pth order
        results[0:-1,8:10] = np.log2(results[0:-1,[3,5]]/results[1:,[3,5]])
        headers = ['$n_{cells}$','$h$','$E_f(h)$','$p_f$','$E_s(h)$','$p_s$']
        results = results[:,[0,1,3,8,5,9]]

        try:
            from tabulate import tabulate
            strt = tabulate(results,headers=headers,tablefmt="latex",floatfmt="4.3e")
            dstfile = os.path.join(self.reportdir,'table_errors.tex')
            if self.debug:
                print strt
        except:
            print 'Unable to create table of errors, have you done pip install tabulate'
            pass

        rate = np.zeros([2,2])
        x = results[np.where(results[:,3]>=self.pth),0]
        r = results[np.where(results[:,3]>=self.pth),2].T
        rate[0,:] = self.linear_fit(np.log10(x),np.log10(r))

        x = results[np.where(results[:,5]>=self.pth),0]
        r = results[np.where(results[:,5]>=self.pth),4].T
        rate[1,:] = self.linear_fit(np.log10(x),np.log10(r))

        try:
            strt = strt + '\n'+ tabulate(rate,headers=['$p$','$c$'],tablefmt="latex",floatfmt="4.3e")
            f = open(dstfile,'w`')
            f.write(strt)
            f.close()
            if self.debug:
                print strt
        except:
            pass

        self.summary['convergence'] = results
        self.summary['log_convergence'] = np.log10(results)
        self.summary['linear_fit'] = rate

        np.save(os.path.join(self.savedir,'convergence_f'+str(self.frame)),results)
        savemat(os.path.join(self.savedir,'summary_f'+str(self.frame)),self.summary)

        # plot results
        if self.plot:
            save_res_name = '_convergence_'
            x = results[:,0]
            for n,case in enumerate(self.convergence_cases):
                y = results[:,2*(n+1)]
                self.plot_convergence(x,y,rate=rate[n,:],log=True,figname=save_res_name,case=case,tolatex=True)