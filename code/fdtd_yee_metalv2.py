# GPL3, Copyright (c) Max Hofheinz, GEGI, UdeS, 2021

import numpy, fiddle, matplotlib, matplotlib.pyplot
from scipy import constants

EFIELD = 0
BFIELD = 1
ENERGY_DENSITY = 2
POYNTING = 3
METAL = 4
NORM=3
DECIBEL=4

def curl_E(E):
    """
    Calculate curl of E
    :param E: E field on Yee grid positions. E is a 4-d array with indices (x, y, z, field_compoent)
    :return: curl of E at Yee grid positions of B field.
    """

    curl_E = numpy.zeros(E.shape)
    curl_E[:, :-1, :, 0] += E[:, 1:, :, 2] - E[:, :-1, :, 2]
    curl_E[:, :, :-1, 0] -= E[:, :, 1:, 1] - E[:, :, :-1, 1]

    curl_E[:, :, :-1, 1] += E[:, :, 1:, 0] - E[:, :, :-1, 0]
    curl_E[:-1, :, :, 1] -= E[1:, :, :, 2] - E[:-1, :, :, 2]

    curl_E[:-1, :, :, 2] += E[1:, :, :, 1] - E[:-1, :, :, 1]
    curl_E[:, :-1, :, 2] -= E[:, 1:, :, 0] - E[:, :-1, :, 0]
    return curl_E

def curl_B(B):
    """
    Calculate curl of B
    :param B: B field on Yee grid positions. B is a 4-d array with indices (x, y, z, field_component)
    :return: curl of B at Yee grid positions of E field.
    """
    curl_B = numpy.zeros(B.shape)

    curl_B[:,1:,:,0] += B[:,1:,:,2] - B[:,:-1,:,2]
    curl_B[:,:,1:,0] -= B[:,:,1:,1] - B[:,:,:-1,1]

    curl_B[:,:,1:,1] += B[:,:,1:,0] - B[:,:,:-1,0]
    curl_B[1:,:,:,1] -= B[1:,:,:,2] - B[:-1,:,:,2]

    curl_B[1:,:,:,2] += B[1:,:,:,1] - B[:-1,:,:,1]
    curl_B[:,1:,:,2] -= B[:,1:,:,0] - B[:,:-1,:,0]
    return curl_B

def poynting(E, B):
    """Calculate Poynting vector from E and B"""
    return numpy.cross(E,B)/constants.mu_0
        



def refine_metal(is_metal):
    
    """
    Convert metallicity from a scalar field to a vector field to make
    sure the metal is terminated with tangential component on the Yee
    grid. 

    """
    # convert to vector field
    is_metal = is_metal[:,:,:,None] + numpy.array([False, False, False])[None,None,None,:]
    # set normal components to 0 at the upper limit of the metal in each direction
    # on the lower limit the condition is automatically verified thanks to the Yee grid
    is_metal[:-1,:,:,0] &= is_metal[1:,:,:,0]
    is_metal[:,:-1,:,1] &= is_metal[:,1:,:,1]
    is_metal[:,:,:-1,2] &= is_metal[:,:,1:,2]
    return numpy.nonzero(is_metal)
    


def timestep(E, B, c, source_pos, source_val, metal_pos):
    """
    Propagate E and B field by 1 full time step
    :param E: renormalized electric field  (4-d array with indices (x, y, z, field_component)) on Yee grid
    :param B: renormalized magnetic field  (4-d array with indices (x, y, z, field_component)) on Yee grid
    :param c: renormalized speed of light in units of space_step/time_step, must be < 1/sqrt(3)
    :param source_pos: positions of source terms
    :param source_val: values of source terms
    :return: renormalized electric field, renormalized magnetic field

    RENORMALIZATION:
    In order to simplify the code we use renormalized values of c, E and B. This avoids having to define
    the time_step and space_step explicitly and include fundamental constants.

    If you want the quantitities in SI units and the have given the current density source_val in A/m^2,
    you have to multiply the fields with the following values to obtain SI units:
    Electric field E in V/m:  E * time_step / epsilon_0
    Magnetic flux density B in T: B / c * time_step**2 / epsilon_0 / space_step = B * time_step * sqrt(mu_0/epsilon_0)
    Magnetic field H in A/m: B / c * time_step**2 / epsilon_0 / mu_0 / space_step = B * c * space_step

    The speed of light c is given in units of space_step/time_step. To get back the speed of light in m/s:
    speed of light in m/s: c * space_step/time_step
    """
    E += c * curl_B(B)

    E[source_pos] += source_val

    E[metal_pos] = 0

    B -= c * curl_E(E)

    return E, B


class WaveEquation:
    """
    Wrapper for live plotting. The __call__ method will be called at regular intervals
    """

    def __init__(self, s, space_step, time_step, c, source, metal, center, radius, int_start,int_stop):
        """
        :param s: 3-tuple giving the shape of the grid
        :param c: renormalized speed of light must be < 1/sqrt(3)
        :param source: function defining the source terms. Takes the time index as input
                       and returns source_pos and source_val (see timestep)
        """
        s = s + (3,)
        self.E = numpy.zeros(s)
        self.B = numpy.zeros(s)
        self.c = c
        self.source = source
        self.index = 0
        self.metal = refine_metal(metal)
        self.space_step = space_step
        self.time_step = time_step
        self.center = center
        self.int_start = int_start
        self.int_stop = int_stop
        self.phi = numpy.linspace(0,2*numpy.pi,200)
        self.cp = numpy.round(numpy.cos(self.phi)*radius).astype(int)
        self.sp = numpy.round(numpy.sin(self.phi)*radius).astype(int)
        self.z = numpy.zeros(self.phi.shape,dtype=int)
        self.pattern_xy = 0
        self.pattern_yz = 0
        self.pattern_xz = 0
        self.npattern = 0


    def update_radiation_pattern(self):
        def pattern(x,y,z):
            return numpy.sum(poynting(self.E[self.center[0]+x, self.center[1]+y, self.center[2]+z],
                                   self.B[self.center[0]+x, self.center[1]+y, self.center[2]+z])*
                             numpy.transpose([x,y,z]),axis=-1)/(x**2+y**2+z**2)
            
        self.pattern_xy = self.pattern_xy + pattern(self.cp,self.sp,self.z)
        self.pattern_yz = self.pattern_yz + pattern(self.z,self.cp,self.sp)
        self.pattern_xz = self.pattern_xz + pattern(self.cp,self.z,self.sp)
        self.npattern += 1

    def plot_radiation_patterns(self):
        def plot_radiation_pattern(pattern, name):
            fig = matplotlib.pyplot.figure()
            ax = fig.add_subplot(projection='polar')
            ax.set_title(name)
            ax.plot(self.phi,10*numpy.log10(pattern/self.npattern))
        
        plot_radiation_pattern(self.pattern_xy,"Plan XY")
        plot_radiation_pattern(self.pattern_yz,"Plan YZ")
        plot_radiation_pattern(self.pattern_xz,"Plan XZ")
        matplotlib.pyplot.show()
        
        
        
        

    def __call__(self, figure, field, component, slice, slice_index, initial=False):
        """
        Perform one time step and plot selected field component
        :param figure: figure object on which to plot
        :param field_component: field component to plot:
        0->Ex, 1->Ey, 2->Ez, 3->Bx 4->By, 5->Bz, 6->Sx, 7->Sy, 8->Sy, 9: Metal
        :param slice: coordinate that will be fixed for 2d plotting 0->x, 1->y, 2->z
        :param slice_index: value of the fixed coordiante
        :param initial: boolean, True if the plot needs to be initialized
        :return:
        """
        #update fields
        
        source_pos, source_index = source(self.index)
        self.E, self.B = timestep(self.E, self.B, self.c, source_pos, source_index, self.metal)
        if self.index >= self.int_start and self.index < self.int_stop:
            self.update_radiation_pattern()
            if self.index+1 >= self.int_stop:
                self.plot_radiation_patterns()
            
        self.index += 1
        lims=1
        if field == EFIELD:
            toplot = self.E
            lims = 1e-2
        elif field == BFIELD:
            toplot = self.B
            lims = 1e-2
        elif field == ENERGY_DENSITY:
            toplot = 0.5*(numpy.sum(self.E**2,axis=-1) + numpy.sum(self.B**2,axis=-1))
            lims = 1e-4
        elif field == POYNTING:
            lims = 5
            toplot = poynting(self.E, self.B)

        elif field == METAL:
            toplot = numpy.zeros(self.E.shape)
            lims = 1
            toplot[self.metal] = 1
            
        if slice == 0:
            toplot = toplot[slice_index, :, :]
            labels = 'yz'
        elif slice == 1:
            toplot = toplot[:, slice_index, :]
            labels = 'xz'
        elif slice == 2:
            toplot = toplot[:, :, slice_index]
            labels = 'xy'
        is_vector_field = len(toplot.shape)==3
        if is_vector_field:
            if component < 3:
                toplot = toplot[:,:,component]
            else:
                # calculate
                is_vector_field = False
                toplot = numpy.sqrt(numpy.sum(toplot**2,axis=-1))
        # imshow expects y coordinate first
        toplot = toplot.transpose()
        if is_vector_field:
            norm = matplotlib.colors.CenteredNorm(halfrange=lims)
            cmap = matplotlib.cm.bwr
        else:
            cmap = matplotlib.cm.cividis
            if component == DECIBEL:
                norm = matplotlib.colors.LogNorm(vmin=lims*1e-3,vmax=10*lims)
            else:
                norm = matplotlib.colors.Normalize(vmin=0,vmax=lims)
        if initial:
            self.axes = figure.add_subplot(111)
            self.image = self.axes.imshow(norm(toplot),norm=None, vmin=0,vmax=1, cmap=cmap,
                                          extent=(-0.5*self.space_step,
                                                  (toplot.shape[0]-0.5)*self.space_step,
                                                  -0.5*self.space_step,
                                                  (toplot.shape[1]-0.5)*self.space_step))
        else:
            self.image.set_data(norm(toplot))
            self.image.set_cmap(cmap)
        self.axes.set_xlabel(labels[0] + ' (m)')
        self.axes.set_ylabel(labels[1] + ' (m)')
        self.axes.set_title('index %d t = %.0f ps' % (self.index, self.time_step * self.index*1e12))
        

def cantenna(grid_dim, base, rmin, height, thickness,wallDistanceFromSource =0 ):
    """
    Draw a cantenna
    :param grid_dim: 3-tuple defining the dimensions of the grid
    :base: position of the bottom of the can (in voxels)
    :rmin: internal radius of the can (in voxels)
    :height: internal height of the can (in voxels)
    :thickness: wall thicknes of the can (in voxels)
    """
    x = numpy.arange(grid_dim[0])
    y = numpy.arange(grid_dim[1])
    z = numpy.arange(grid_dim[2])
    r = numpy.sqrt((x[:,None]-base[0])**2+(y[None,:]-base[1])**2)
    cylinder = ((r >= rmin) & (r <= (rmin+thickness)))[:,:,None] & ((z <= base[2] + height) & (z >= base[2]))[None,None,:]
    bottom = (r <= (rmin+thickness))[:,:,None] & ((z <= base[2]) & (z >= base[2]-thickness))[None,None,:]
    wall = (r <= (rmin+thickness))[:,:,None] & ((z <= base[2] + wallDistanceFromSource) & (z >= base[2]+wallDistanceFromSource-thickness))[None,None,:]
    return cylinder | bottom | wall


    
    


if __name__ == "__main__":
    put_cantenna=False # put cantenna around dipole antenna or not
    n = 100 # grid size n x n x n
    f = 2.4e9 # source frequency
    time_step = 1.0/f/200 # time step in s
    c = 0.1 # renormalized speed of light in voxel/iteration, must be < 1/sqrt(3)
    space_step = constants.c * time_step / c
    cantenna_radius = 0.095/2 / space_step
    cantenna_height = 0.133 / space_step
    antenna_from_bottom = 0.01 / space_step
    cantenna_thickness = 3# make walls at least 3 voxels thick, otherwise there might be discretization problems making the can leaky
    cantenna_bottom = n//2 - cantenna_height
    if put_cantenna:
        radiation_diagram_start = 400
        source_pos = int(cantenna_bottom + antenna_from_bottom)
    else:
        radiation_diagram_start = 300
        source_pos = n//2
    radiation_diagram_center = [n//2,n//2,n//2]
    radiation_diagram_radius = 0.35*n

    radiation_diagram_stop = radiation_diagram_start + 1./f/time_step
    

    def source(index):
        """
        source term of electric field (current)
        :param index: time step number
        :return: source_pos, source_val
        source_pos: coordinates of sources
        source_val: corresponding current values
        """
        #return current source in x direction (last index=0) at coordinates (50,50,20)
        return ([50], [50], [source_pos], [0]), 0.1*numpy.sin(2*numpy.pi*f*time_step * index)
    
    # size of the space grid
    dims = (n,n,n)
    
    if put_cantenna:
        # simulation with cantenna
        w = WaveEquation(dims, space_step, time_step, c, source,
                         cantenna(dims, (n//2,n//2,cantenna_bottom),cantenna_radius,
                                  cantenna_height,cantenna_thickness),
                         radiation_diagram_center,radiation_diagram_radius,
                         radiation_diagram_start,radiation_diagram_stop,30)
    else:
        # simulation without metal objects
        w = WaveEquation(dims, space_step, time_step, c, source, numpy.zeros(dims,dtype=bool),
                         radiation_diagram_center, radiation_diagram_radius,
                         radiation_diagram_start,radiation_diagram_stop)

    fiddle.fiddle(w, [('field',{'E':EFIELD,'B':BFIELD,'Energy density':ENERGY_DENSITY, 'Poynting':POYNTING, 'Metal':METAL},'E'),
                       ('component',{'X':0, 'Y':1, 'Z':2,'norm':NORM,'dB':DECIBEL},'norm'),
                      ('slice',{'XY':2,'YZ':0,'XZ':1},'XZ'),
                      ('slice index',0,n-1,50,1)], update_interval=0.01)


