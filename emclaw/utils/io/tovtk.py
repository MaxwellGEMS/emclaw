from inout import IO
import os
import numpy as np
import scipy as sp
import scipy.special as sps
from scipy.integrate import simps as sq
import numpy.linalg as npl

def QuadField1D(solution,quad):
    #diff = npl.norm(solution.q[1]-q_old[1],1)
    quad[0,solution.frame] = solution.t
    quad[1,solution.frame] = sq(solution.q[0]*solution.q[1],solution.x.centers)
    quad[2,solution.frame] = sq(np.sqrt(solution.q[0]**2 + solution.q[1]**2),solution.x.centers)    
    return quad

def QuadField2D(solution,quad):
    #diff = npl.norm(solution.q[1]-q_old[1],1)
    quad[0,solution.frame] = solution.t
    quad[1,solution.frame] = sq(np.sqrt(solution.q[1]*solution.q[2]**2 + (-solution.q[0]*solution.q[2])**2),solution.x.centers)
    quad[2,solution.frame] = sq(np.sqrt(solution.q[0]**2 + solution.q[1]**2 + solution.q[2]**2),solution.x.centers)    
    return quad

def Poyinting1D(solution):
    S = np.zeros([1,len(solution.x.centers)])
    S = solution.q[0]*solution.q[1]
    return S

def Poyinting2D(solution):
    S = np.zeros([2,len(solution.x.centers),len(solution.y.centers)])
    S[0,:,:] = solution.q[1]*solution.q[2]
    S[1,:,:] = -solution.q[0]*solution.q[2]


if __name__ == "__main__":
    import sys
    path = sys.argv[1]
    num_frames = int(sys.argv[2])
    print 'going to path:', path
    print 'number of frames:', num_frames
    quad = np.zeros([3,num_frames+1])
    file_name = os.join.path(path,'quad.txt')
    for i in range(0,num_frames+1):
        print i
        sol = IO()
        sol.path = path
        sol.frame = i
        sol.read_petsc()
        sol.q_to_vtk()
        postcalc(sol,quad)
        q_old = sol.q.copy()

    np.savetxt(file_name,quad)
    os.remove('petclaw.log')
    os.remove('pyclaw.log')
    os.remove('inout.pyc')