#import base
import sys, os, csv, string
from pylab import *
from pprint import *

#import Qt 
from PySide.QtCore import Signal
from PySide.QtCore import Slot
from PySide.QtCore import *
from PySide.QtGui import *
from PySide import QtCore, QtGui

#spectra fitting toolkit
from spectra_fitting import *

class ParameterSlider(QAbstractSlider):

  changed = Signal()

  def __init__(self, var, parent=None):
    super(ParameterSlider,self).__init__(parent)
    self.parent = parent

    self.ignore_signals = False

    self.varLabel = QLabel(var)

    self.value = QLineEdit()
    self.value.textChanged.connect(self.updateSlider)
    self.value.textChanged.connect(self.changed)

    self.slider = QSlider(Qt.Horizontal)
    self.slider.setMinimum(0)
    self.slider.setMaximum(1000)
    self.slider.setSingleStep(1)
    self.slider.valueChanged[int].connect(self.sliderChanged)

    inLabel = QLabel('\in (')
    commaLabel = QLabel(',')
    endLabel = QLabel(')')

    self.lower = QLineEdit('0.0')
    self.lower.textChanged.connect(self.updateSlider)
    self.lower.textChanged.connect(self.changed)
    self.upper = QLineEdit('1.0')
    self.upper.textChanged.connect(self.updateSlider)
    self.upper.textChanged.connect(self.changed)

    self.maxValue = 1000.0
    self.minValue = 0.0

    upperBox = QHBoxLayout()
    mainLayout = QVBoxLayout()
    upperBox.addWidget(self.varLabel)
    upperBox.addWidget(self.value)
    upperBox.addWidget(inLabel)
    upperBox.addWidget(self.lower)
    upperBox.addWidget(commaLabel)
    upperBox.addWidget(self.upper)
    upperBox.addWidget(endLabel)

    mainLayout.addLayout(upperBox)
    mainLayout.addWidget(self.slider)
    self.setLayout(mainLayout)

  def sliderChanged(self, value):
    if self.ignore_signals:
      return

    self.ignore_signals = True
    try:
      self.value.setText(str(self.minValue + value * (self.maxValue-self.minValue)/1000) )
    except:
      pass
    self.ignore_signals = False

  def updateSlider(self):
    if self.ignore_signals:
      return

    self.ignore_signals = True
    try:
      self.minValue = float(self.lower.text())
      self.maxValue = float(self.upper.text())
      self.slider.setValue( int( (float(self.value.text())-self.minValue) * 1000/(self.maxValue-self.minValue)))
    except:
      pass
    self.ignore_signals = False

  def setMinimum(self, value):
    self.minValue = value
    self.lower.setText(str(value))

  def setMaximum(self, value):
    self.maxValue = value
    self.upper.setText(str(value))

  def getValue(self):
    return float(self.value.text())

  def setValue(self, value):
    self.value.setText(str(value))

  def getRange(self):
    return r_[float(self.lower.text()), float(self.upper.text())]

class ParameterDialog(QDialog):
    def __init__(self, fit_object, parent=None):
      super(ParameterDialog,self).__init__(parent)

      self.parent = parent
    
      self.fit_object = fit_object

      spec = fit_object.get_spec()
      self.parametersGroup = QGroupBox("Edit Parameters")
 
      self.paramSliders = []

      parametersLayout = QVBoxLayout()
      for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
          paramSlider = ParameterSlider(var)
          self.paramSliders.append(paramSlider)
          paramSlider.setValue(val)
          paramSlider.setMinimum(range[0])
          paramSlider.setMaximum(range[1])
          paramSlider.changed.connect(self.update)
          #paramSlider.changed.connect(fit_object.optimization_window.update)

          parametersLayout.addWidget(paramSlider)

      self.parametersGroup.setLayout(parametersLayout)
      mainLayout = QVBoxLayout()
      mainLayout.addWidget(self.parametersGroup)
      self.setLayout(mainLayout)
      self.setWindowTitle(fit_object.name)

    def update(self):
      spec = self.fit_object.get_spec()
      values = []
      ranges = []

      for paramSlider in self.paramSliders:
        values.append(paramSlider.getValue())
        ranges.append(paramSlider.getRange())

      values = array(values)

      spec['values']=values
      spec['ranges']=ranges

      self.fit_object.set_spec(spec)
      self.parent.on_show()

