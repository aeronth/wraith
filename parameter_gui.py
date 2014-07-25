#import base
import sys, os, csv, string
from pylab import *
#from pprint import *

#import Qt
#from PySide.QtCore import Signal
#from PySide.QtCore import Slot
from PySide.QtCore import *
from PySide.QtGui import *
#from PySide import QtCore, QtGui

#spectra fitting toolkit
from spectra_fitting import *

class ParameterSlider(QAbstractSlider):
    """ dynamically generates controls to tweak float parameters
        in dialogs with 'n' adjustable parameters 
    """
    changed = Signal()

    def __init__(self, var, parent=None):
        super(ParameterSlider,self).__init__(parent)
        self.parent = parent

        self.ignore_signals = False

        self.varLabel = QLabel(var)

        #replace QLineEdit() with spin box DFO
        self.minSlider = 0
        self.rangeSlider = 1000
        self.inc = 0.01 #for spin boxes
        self.range_inc = 0.1
        self.maxValue = self.minSlider + self.rangeSlider
        self.minValue = self.minSlider

        #parameter value and limits
        self.value = QDoubleSpinBox()
        self.value.setAccelerated(True )
        self.value.setCorrectionMode(QAbstractSpinBox.CorrectToNearestValue)
        self.value.setKeyboardTracking(False)
        self.value.setRange(self.minValue,self.maxValue)
        self.value.valueChanged.connect(self.setSlider)
        self.value.valueChanged.connect(self.changed)

        self.lower = QDoubleSpinBox()
        self.lower.setKeyboardTracking(False)
        self.lower.setValue(self.minValue)
        self.lower.valueChanged.connect(self.setSlider)
        self.lower.valueChanged.connect(self.changed)
        self.upper = QDoubleSpinBox()
        self.upper.setKeyboardTracking(False)
        self.upper.setValue(self.maxValue)     
        self.upper.valueChanged.connect(self.setSlider)
        self.upper.valueChanged.connect(self.changed)
        
        self.setScrollStep()
        #parameter value slider and decorations
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setTracking(False)  #update after slider released
        self.slider.setRange(self.minSlider, self.minSlider + self.rangeSlider)
        self.slider.setSingleStep(int(self.inc * self.rangeSlider))
        self.slider.valueChanged.connect(self.sliderChanged)

        #per-parameter layout
        inLabel = QLabel('in (')
        commaLabel = QLabel(',')
        endLabel = QLabel(')')

        upperBox = QHBoxLayout()
        upperBox.addWidget(self.varLabel)
        upperBox.addWidget(self.value)
        upperBox.addWidget(inLabel)
        upperBox.addWidget(self.lower)
        upperBox.addWidget(commaLabel)
        upperBox.addWidget(self.upper)
        upperBox.addWidget(endLabel)

        mainLayout = QVBoxLayout()
        mainLayout.addLayout(upperBox)
        mainLayout.addWidget(self.slider)
        self.setLayout(mainLayout)

    def sliderChanged(self, valueSlider):
        if self.ignore_signals:
            return

        self.ignore_signals = True
        try:
            #self.value.setText(str(self.minValue + value * (self.maxValue-self.minValue)/1000) )
            self.value.setValue(self.minValue + valueSlider * (self.maxValue-self.minValue)/self.rangeSlider )
        except:
            pass
        self.ignore_signals = False

    def updateSlider(self):
        if self.ignore_signals:
            return

        self.ignore_signals = True
        try:
            #self.minValue = self.lower.value()
            #self.maxValue = self.upper.value()
            #self.slider.setValue( int( (float(self.value.text())-self.minValue) * 1000/(self.maxValue-self.minValue)))
            #self.slider.setValue( int( (self.value.getValue()-self.minValue) * self.rangeSlider/(self.maxValue-self.minValue)))
            self.setSlider()
        except:
            pass
        self.ignore_signals = False

    def setMinimum(self, value):
        self.minValue = value
        self.lower.setValue(value)
        self.value.setRange(self.minValue,self.maxValue)
        self.setScrollStep()
        self.setSlider()

    def setMaximum(self, value):
        self.maxValue = value
        self.upper.setValue(value)
        self.value.setRange(self.minValue,self.maxValue)
        self.setScrollStep()
        self.setSlider()
        
    def setScrollStep(self):
        step = self.inc*abs(self.maxValue - self.minValue)
        self.value.setSingleStep(step)
        if step < 1:
            decimal = ceil(-log10(step))
            self.value.setDecimals(decimal)
        step = self.range_inc*abs(self.maxValue - self.minValue)
        self.lower.setSingleStep(step)
        self.upper.setSingleStep(step)
        if step < 1:
            decimal = ceil(-log10(step))
            self.lower.setDecimals(decimal)
            self.upper.setDecimals(decimal)      

    def getValue(self):
        return self.value.value()

    def setValue(self, value):
        self.value.setValue(value)
        self.setSlider()
        
    def setSlider(self):
        self.slider.setValue( int( self.minSlider + (self.getValue()-self.minValue) * self.rangeSlider/(self.maxValue-self.minValue)))

    def getRange(self):
        return r_[self.lower.value(), self.upper.value()]

class ParameterDialog(QDialog):
    """ Dymanically creates dialog with ParameterSliders for each float argument
    """
    def __init__(self, fit_object, parent=None):
        super(ParameterDialog,self).__init__(parent)

        self.parent = parent

        self.fit_object = fit_object

        spec = fit_object.get_spec()
        self.parametersGroup = QGroupBox("Edit Fit Parameters")

        self.paramSliders = []

        parametersLayout = QVBoxLayout()
        for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
            paramSlider = ParameterSlider(var)
            
            paramSlider.setMinimum(range[0])
            paramSlider.setMaximum(range[1])
            paramSlider.setValue(val)            
           
            self.paramSliders.append(paramSlider)
            #paramSlider.changed.connect(fit_object.optimization_window.update)
            parametersLayout.addWidget(paramSlider)
        
        for paramSlider in self.paramSliders:
            paramSlider.changed.connect(self.update)
        
        #Maybe?
        #self.update()
        
        self.parametersGroup.setLayout(parametersLayout)
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.parametersGroup)
        self.setLayout(mainLayout)
        self.setWindowTitle(fit_object.name)    #put treename in spec directory

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
