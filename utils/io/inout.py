# petclaw to vtk

import os
import numpy as np
from petsc4py import PETSc
import pickle
import glob
import shutil

def post_calculation():
        pass

class IO(object):

    def read_petsc(self):

        if hasattr(self, 'frame'): frame = self.frame
        if hasattr(self,'file_prefix'): file_prefix = self.file_prefix
        if hasattr(self, 'path'): path = self.path
        if hasattr(self, 'write_aux'): write_aux = self.write_aux
        if hasattr(self, 'write_aux'): read_aux = self.read_aux
        if hasattr(self, 'write_p'): write_p = self.write_p

        pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
        viewer_filename = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)
        aux_viewer_filename1 = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)
        aux_viewer_filename2 = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(0).zfill(4)
        if os.path.exists(aux_viewer_filename1):
             aux_viewer_filename = aux_viewer_filename1
        else:
             aux_viewer_filename = aux_viewer_filename2

        pickle_file = open(pickle_filename,'rb')

        # this dictionary is mostly holding debugging information, only nstates is needed
        # most of this information is explicitly saved in the individual patches
        value_dict = pickle.load(pickle_file)
        nstates = value_dict['nstates']
        num_dim = value_dict['num_dim']
        num_aux = value_dict['num_aux']
        num_eqn = value_dict['num_eqn']

        self.__setattr__('num_dim',num_dim)
        self.__setattr__('num_aux',num_aux)
        self.__setattr__('num_eqn',num_eqn)
        # now set up the PETSc viewer (assuming binary)
        viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)

        patches = []
        for m in xrange(nstates):
            patch_dict = pickle.load(pickle_file)

            level = patch_dict['level']
            names = patch_dict['names']
            lower = patch_dict['lower']
            n = patch_dict['num_cells']
            d = patch_dict['delta']
            from clawpack import petclaw
            dimensions = []
            for i in xrange(num_dim):
                dimensions.append(
                    petclaw.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
            patch = petclaw.Patch(dimensions)
            self.__setattr__('_patch',patch)

            if num_dim==1:
                self.__setattr__('x',patch.x)
            elif num_dim==2:
                self.__setattr__('x',patch.x)
                self.__setattr__('y',patch.y)
            elif num_dim == 3:
                self.__setattr__('y',patch.y)
                self.__setattr__('z',path.z)


            self.__setattr__('num_cells',patch.num_cells_global)
            claw = petclaw.State(patch,num_eqn,num_aux) ##
            self.__setattr__('_claw',claw)
            self.t = value_dict['t']
            self.problem_data = value_dict['problem_data']
            self.nstates = value_dict['nstates']
            self._claw.gqVec.load(viewer)
            if read_aux:
                self._claw.gauxVec.load(aux_viewer)
                self.__setattr__('aux',self._claw.aux)
            
            self.__setattr__('q',self._claw.q)
            self.__setattr__('frame',frame)
            self.__setattr__('file_prefix',file_prefix)
            self.__setattr__('read_aux', read_aux)
            self.__setattr__('write_aux', write_aux)
            self.__setattr__('write_p', write_p)
            self.__setattr__('path', path)

        return self

    def write_vtk(self):
        
        if hasattr(self, 'frame'): frame = self.frame
        if hasattr(self,'file_prefix'): file_prefix = self.file_prefix
        if hasattr(self, 'path'): path = self.path
        if hasattr(self, 'write_aux'): write_aux = self.write_aux
        if hasattr(self, 'write_p'): write_p = self.write_p
        
        if hasattr(self, 'q'):
            viewer_filename = os.path.join(path, file_prefix+str(frame).zfill(4)+'.vtk')
            viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE, format=PETSc.Viewer.Format.ASCII_VTK)

            self.gqVec.view(viewer)
        
            if write_aux:
                self.gauxVec.view(aux_viewer)
        
            viewer.flush()
            viewer.destroy()
            if write_aux:
                aux_viewer.flush()
                aux_viewer.destroy()

    def q_to_vtk(self):

        nx,ny = self.num_cells

        coordinates = [self.x.centers,
               self.y.centers,
               np.ones(1),
               ]
        
        dimensions = (nx, ny, 1) 
        
        if hasattr(self, 'frame'): frame = self.frame
        if hasattr(self,'file_prefix'): file_prefix = self.file_prefix
        if hasattr(self, 'path'): path = self.path
        
        scalars = [("Q1", self.q[0]),
           ("Q2", self.q[1]),
           ("Q3", self.q[2])]

        vectors = []
        title = 'VTK Data'
        
        filename = os.path.join(path, file_prefix+str(frame).zfill(4)+'.vtk')
        fh = open(filename, 'wb')
        fh_write = lambda s: fh.write(s.encode('ascii'))

        header = '# vtk DataFile Version %d.%d'
        version = (2, 0)
        fh_write(header % version)
        fh_write('\n')
        title = title
        fh_write(title[:255])
        fh_write('\n')

        format = 'BINARY'
        fh_write(format)
        fh_write('\n')

        dataset_type = 'RECTILINEAR_GRID'
        fh_write('DATASET %s' % dataset_type);
        fh_write('\n')
        fh_write('DIMENSIONS %d %d %d' % dimensions)
        fh_write('\n')
        for X, array in zip("XYZ", coordinates):
            label = X+'_COORDINATES'
            fh_write('%s %s %s' % (label, len(array), 'double'))
            fh_write('\n')
            array.astype('>d').tofile(fh)
            fh_write('\n')

        data_type = 'POINT_DATA'
        fh_write('%s %d' % (data_type, np.prod(dimensions)))
        fh_write('\n')

        for i, (name, array) in enumerate(scalars):
            attr_type = 'SCALARS'
            attr_name = name or (attr_type.lower() + str(i))
            attr_name = attr_name.replace(' ', '_')
            fh_write('%s %s %s' %(attr_type, attr_name, 'double'))
            fh_write('\n')
            lookup_table = 'default'
            lookup_table = lookup_table.replace(' ', '_')
            fh_write('LOOKUP_TABLE %s' % lookup_table)
            fh_write('\n')
            array.astype('>d').tofile(fh)
            fh_write('\n')

        for i, (name, array) in enumerate(vectors):
            attr_type = 'VECTORS'
            attr_name = name or (attr_type.lower() + str(i))
            attr_name = attr_name.replace(' ', '_')
            fh_write('%s %s %s' %(attr_type, attr_name, 'double'))
            fh_write('\n')
            array.astype('>d').tofile(fh)
            fh_write('\n')

        fh.flush()
        fh.close()

    def copy_files(self,src_file,dst_file):
        pass

    def __init__(self,frame=0,file_prefix='claw',path='./',write_aux=False,write_p=False,read_aux=False):

        self.frame=frame
        self.file_prefix=file_prefix 
        self.path=path
        self.read_aux = read_aux
        self.write_aux = write_aux
        self.write_p = write_p
        self.write_postprocess = False
        self.postprocess = post_calculation()

