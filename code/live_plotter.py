# Copyright Max Hofheinz, Universite de Sherbrooke. max.hofheinz@usherbrooke.ca
# License: GPL2
# Date 2018/07/18


import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

class ParamFiddle(tk.Frame):
    """
    A widget controlling how a plot parameter is handled. It allows to
    chose in to be plotted along X, Y or to have a fixed value.  The
    callback routine provided in the constructor is called whenever
    the state changes
    """
    
    def onValueChange(self, adj):
        newvalue = self.scale.get()
        if self.lastvalue != newvalue:
            self.lastvalue = newvalue
            self.callback(self.paramIndex, newvalue)
    
    def getValue(self):
        return self.lastvalue


    def __init__(self, parent, callback, paramIndex,
                 paramName, vMin, vMax, vStart=None):
        self.parent = parent
        self.vMin = vMin
        self.vMax = vMax
        if vStart is None:
            self.lastvalue = 0.5*(vMin+vMax)
        else:
            self.lastvalue = vStart
        self.paramIndex = paramIndex
        self.callback = callback

        tk.Frame.__init__(self, parent)
        self.state = tk.StringVar()
        self.label = tk.Label(self, text=paramName,width=18)
        self.label.pack(side=tk.LEFT)
        self.scale = tk.Scale(self, from_ = vMin, to=vMax,
                                   resolution = -1, orient = tk.HORIZONTAL,
                                   showvalue=1, command = self.onValueChange)
        self.scale.set(self.lastvalue)
        self.scale.pack(side=tk.LEFT,expand=1,fill=tk.X)


class FiddlePlotter(tk.Frame):
    def __init__(self, parent, parameters, update_interval=False):
        tk.Frame.__init__(self, parent)
        self.parameters = parameters
        #color = tk.Style().lookup('TFrame','background')
        color = None
        self.figure = Figure(figsize=(6,4), dpi=60, linewidth=0,facecolor=color)
        self.canvas = FigureCanvasTkAgg(self.figure,master=self)
     #   canvas.show()
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH,expand=1)
        self.paramControls = [ParamFiddle(self, self.onParamControl, i, *p) \
                              for i,p in enumerate(parameters)]
        for p in self.paramControls:
            p.pack(side=tk.TOP, fill=tk.X)
        self.paramvalues = [p.getValue() for p in self.paramControls]    
        toolbar = NavigationToolbar2Tk(self.canvas, self)
        toolbar.update()
        toolbar.pack(side=tk.TOP, fill=tk.X)
        self.updateID = None
        self.update_interval = False
        if self.initial() and update_interval:
            self.update_interval = int(update_interval*1e3)
            self.updateID = self.after(self.update_interval, self._update)

    def onParamControl(self, paramIndex, value):
        if self.updateID is not None:
            self.after_cancel(self.updateID)
        self.paramvalues[paramIndex] = value
        self._update()
        
    def update(self):
        """
        update plot after parameter change or preiodically
        """
        pass
        
    def _update(self):
        self.update()
        self.canvas.draw()
        if self.update_interval:
            self.updateID = self.after(self.update_interval, self._update)


     

    def initial(self):
        pass


    def final(self):
        pass

    def destroy(self):
        if self.updateID is not None:
            self.after_cancel(self.updateID)
        try:
            self.final()
        finally:
            tk.Frame.destroy(self)


def open_plotter(PlotClass, parameters, w=800,h=600,title="Fiddle Plotter", update_interval=False):
    window = tk.Tk()
    window.geometry("%dx%d" % (w,h))
    window.title(title)
    plotter = PlotClass(window, parameters, update_interval=update_interval)
    plotter.pack(expand=tk.YES,fill=tk.BOTH)
    window.mainloop()
    

    
class ListBoxChoice(object):
    def __init__(self, master=None, title=None, message=None, list=[]):
        self.master = master
        self.value = None
        self.list = list[:]

        self.modalPane = tk.Toplevel(self.master)

        self.modalPane.transient(self.master)
        self.modalPane.grab_set()

        self.modalPane.bind("<Return>", self._choose)
        self.modalPane.bind("<Escape>", self._cancel)

        if title:
            self.modalPane.title(title)

        if message:
            tk.Label(self.modalPane, text=message).pack(padx=1, pady=1)

        listFrame = tk.Frame(self.modalPane)
        listFrame.pack(side=tk.TOP, padx=5, pady=5,fill=tk.BOTH,expand=True)

        scrollBar = tk.Scrollbar(listFrame)
        scrollBar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listBox = tk.Listbox(listFrame, selectmode=tk.SINGLE)
        self.listBox.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollBar.config(command=self.listBox.yview)
        self.listBox.config(yscrollcommand=scrollBar.set)
        for item in self.list:
            self.listBox.insert(tk.END, item)

        buttonFrame = tk.Frame(self.modalPane)
        buttonFrame.pack(side=tk.BOTTOM,fill=tk.X)

        chooseButton = tk.Button(buttonFrame, text="Choose", command=self._choose)
        chooseButton.pack(side=tk.LEFT)

        cancelButton = tk.Button(buttonFrame, text="Cancel", command=self._cancel)
        cancelButton.pack(side=tk.RIGHT)

    def _choose(self, event=None):
        sel = self.listBox.curselection()
        if len(sel) ==1: 
            self.value = sel[0]
            self.modalPane.destroy()

    def _cancel(self, event=None):
        self.value = None
        self.modalPane.destroy()

    def returnValue(self):
        self.master.wait_window(self.modalPane)
        return self.value



    
