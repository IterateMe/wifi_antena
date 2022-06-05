#!/usr/bin/python3

# Copyright Max Hofheinz, Universite de Sherbrooke: max.hofheinz@usherbrooke.ca
# License: GPLv3
# Date 2018/09/26
# - correct dephasing between E and B for lossy media
# Date 2018/07/14

import live_plotter, numpy, time

class WavePlotter(live_plotter.FiddlePlotter):

    def initial(self):
        """
        Initialize plot (no wave yet)
        """
        self.axes = self.figure.add_axes((0,0,1,1),aspect='equal',frameon=False)
        self.axes.get_xaxis().set_visible(False)
        self.axes.get_yaxis().set_visible(False)
        
        self.t0 = time.time()
        self.omega = 1.0
        self.n = 101
        self.zscale = 0.2
        self.zmax = 10
        self.z = numpy.linspace(0,self.zmax,self.n)
        self.ivec = numpy.arange(0,self.n,5)
        self.axes.plot([0,self.zmax*self.zscale],[0,self.zmax*self.zscale],'k-',lw=3)
        self.E0 = self.axes.plot([0*self.ivec]*2,[0*self.ivec]*2,'r-',lw=3)
        self.H0 = self.axes.plot([0*self.ivec]*2,[0*self.ivec]*2,'b-',lw=3)

        self.Et, = self.axes.plot(0*self.z,0*self.z,'r-',label='E')
        self.Ht, = self.axes.plot(0*self.z,0*self.z,'b-',label='H')
        self.axes.set_xlim(-1,3)
        self.axes.set_ylim(-1,3)
        self.axes.legend(fontsize=20)
        return True
  

    def map(self, x, y, z):
        return y + 0.5 * z, x + 0.5 * z
    
    def update(self):
        """
        Calculates the plane wave and updates the graph periodically and after each parameter change
        """
        #attibute parameters 
        Ex = self.paramvalues[0]
        Ey = self.paramvalues[1] * numpy.exp(self.paramvalues[2] * 2j * numpy.pi / 360)
        gamma =  self.paramvalues[3] + 1.0j* self.paramvalues[4]
        t = time.time() - self.t0
        #calculate elctric field
        phasor = numpy.exp(1j * self.omega * t - gamma * self.z)
        Ex = Ex * phasor
        Ey = Ey * phasor
        # vec_H = 1/eta * vec_ap x vec_E
        #E and B have arbitrary scale so we set mu and omega to 1, then 1j/eta=gamma
        Hy = Ex * gamma / 1.0j
        Hx = -Ey * gamma / 1.0j
        zi= self.zscale*self.z
        for j,i in enumerate(self.ivec):
            self.E0[j].set_data([zi[i], Ey[i].real + zi[i]], [zi[i], Ex[i].real + zi[i]])
            self.H0[j].set_data([zi[i], Hy[i].real + zi[i]], [zi[i], Hx[i].real + zi[i]])
        self.Et.set_data(Ey.real + zi, Ex.real + zi)
        self.Ht.set_data(Hy.real + zi, Hx.real + zi)
        
              

def wave_plot():
    # define parameter sliders: (name, min, max, initial_value)
    parameters = [("Ex",-1,1,1),
                  ("Ey",-1,1,0),
                  ("phase y vs x",-180,180,0),
                  ("alpha",0,1,0),
                  ("beta",0,10,1)]
    live_plotter.open_plotter(WavePlotter, parameters, title="Plane Wave Animation",update_interval=0.05,w=600,h=800)


if __name__ == "__main__":
    wave_plot()
    
