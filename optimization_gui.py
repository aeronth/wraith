#import base
import sys, os, csv, string
from pprint import *

#import Qt 
from PySide.QtCore import Signal
from PySide.QtCore import Slot
from PySide.QtCore import *
from PySide.QtGui import *
from PySide import QtCore, QtGui

#import pylab, which includes scipy and numpy
from pylab import *

#### import project specific items ####

#spectra fitting toolkit
from spectra_fitting import *

class OptimizationWindow(QMainWindow):
  def __init__(self, fit_object, parent=None):
    super(OptimizationWindow,self).__init__(parent)

    self.parent = parent

    self.fit_object = fit_object

    spec = fit_object.get_spec()

    self.main_frame = QWidget()

    mainLayout = QVBoxLayout()
    self.setWindowTitle('Optimization - ' + fit_object.name)

    self.dpi = 100
    self.fig = Figure((6.0, 4.0), dpi=self.dpi, facecolor='w', edgecolor='k')
    self.canvas = FigureCanvas(self.fig)
    self.canvas.setParent(self.main_frame)

    N = size(spec['variables'])+1
    i = 1
    self.axes = {}
    for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
        self.axes[var] = Subplot(self.fig,ceil(N/2.0),2,i)
        i += 1
        self.fig.add_subplot(self.axes[var])
        self.axes[var].set_ylabel("$"+var+"$")
        self.axes[var].set_ylim(range)

    self.opt_axes = Subplot(self.fig,ceil(N/2.0),2,N)
    self.fig.add_subplot(self.opt_axes)

    self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

    mainLayout.addWidget(self.canvas)
    mainLayout.addWidget(self.mpl_toolbar)

    self.main_frame.setLayout(mainLayout)

    self.setCentralWidget(self.main_frame)
    self.show()
    self.max_r = max(self.fit_object.spectrum.data)*1e2

  def update(self):
    spec = self.fit_object.get_spec()
    i=0
    window_size = 2500
    l = len(self.fit_object.optimization_history[i,:])
    for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
        self.axes[var].cla()
        self.axes[var].plot(self.fit_object.optimization_history[i,:])
        i+=1
        self.axes[var].set_ylim(range)
        self.axes[var].set_ylabel("$"+var+"$")
        self.axes[var].set_xlim([l-window_size,l])

    points = self.fit_object.optimization_history[i,:]
    points[points>self.max_r] = self.max_r
    #points[points<0] = 0
    self.opt_axes.cla()
    self.opt_axes.plot(points)
    self.opt_axes.set_ylabel("$\sum R^2$")
    self.opt_axes.set_xlim([l-window_size,l])
    self.opt_axes.set_ylim(mean(points[-window_size:])+r_[-6,6]*std(points[-window_size:]))

    self.canvas.draw()
    self.show()



