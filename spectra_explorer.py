#! /usr/bin/python

import sys, os, csv, string
#from PySide.QtCore import *
from PySide.QtCore import Signal
from PySide.QtCore import Slot
from PySide.QtCore import *
from PySide.QtGui import *
from PySide import QtCore, QtGui

import matplotlib
matplotlib.use('Qt4Agg')
matplotlib.rcParams['backend.qt4']= "PySide"

from pylab import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure
from spectra_fitting import *
from VAMAS import *
from pprint import *
from mpl_toolkits.axisartist import Subplot


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
          paramSlider.changed.connect(fit_object.opt_window.)

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

    N = size(spec['variables'])
    i = 1
    for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
        self.axes[var] = Subplot(self.fig,N,1,i)
        i += 1
        self.fig.add_subplot(self.axes[var])
        self.axes[var].set_title(var)
        self.axes[var].set_ylim(range)
    
    self.opt_axes = Subplot(self.fig,N,1,i)

    self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
    
    mainLayout.addWidget(self.canvas)
    mainLayout.addWidget(self.mpl_toolbar)

    self.main_frame.setLayout(mainLayout)

    self.setCentralWidget(self.main_frame)
    self.show()

  def update(self):
    spec = self.fit_object.get_spec()
    i=0
    for var, val, range in zip(spec['variables', spec['values'], spec['ranges']):
        self.axes[var].cla()
        self.axes[var].plot(optimization_history[:,i])
        i+=1
        self.axes[var].set_ylim(range)

    points = optimization_history[:,i]
    points = points[points<max(spectrum.data)*1e5]
    points = points[points>0]
    self.opt_axes.cla()
    self.opt_axes.plot(optimization_history[:,i])

    self.canvas.draw()
    self.show()


class AverageWindow(QMainWindow):
    def __init__(self, spectrum, parent=None):
      super(AverageWindow,self).__init__(parent)
  
      self.parent = parent
    
      self.spectrum = spectrum
  
      self.main_frame = QWidget()
  
      mainLayout = QVBoxLayout()
      self.setWindowTitle('Average')

      self.dpi = 100
      self.fig = Figure((6.0, 4.0), dpi=self.dpi, facecolor='w', edgecolor='k')
      self.canvas = FigureCanvas(self.fig)
      self.canvas.setParent(self.main_frame)

      self.axes = Subplot(self.fig,1,1,1)
      self.fig.add_subplot(self.axes)
      self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
      
      mainLayout.addWidget(self.canvas)
      mainLayout.addWidget(self.mpl_toolbar)

      self.main_frame.setLayout(mainLayout)

      self.setCentralWidget(self.main_frame)

      spectrum.plot_full_summary(scale=1,axes=self.axes)

    def update(self):
      pass

class Form(QMainWindow):
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        self.setWindowTitle('Interactive XPS Explorer')

        self.files = {}
        #self.data = DataHolder()
        self.series_list_model = QStandardItemModel()
        self.series_list_root = self.series_list_model.invisibleRootItem()
        self.series_list_model.itemChanged.connect(self.update_file_checks)

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        
        self.update_ui()
        self.on_show()

    def filterSelected(self, filter):
        self.filefilter = filter

    def load_file(self, filename=None):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFiles)
        dialog.setNameFilter('VAMAS (*.vms);; BL7.0.1.1 XAS (*.txt);; SSRL (*.dat);; SUPER (*);; All Files (*.*)')
        self.filefilter = 'VAMAS (*.vms)'
        dialog.filterSelected.connect(self.filterSelected)

        if dialog.exec_():
          filenames = dialog.selectedFiles()

        for filename in filenames:
          if filename:
            file = os.path.basename(str(filename))
            self.files[file] = DataHolder()
            self.files[file].load_from_file(str(filename), str(self.filefilter))
            self.fill_series_list()
            self.update_ui()

    def fill_series_list(self):
        self.series_list_model.clear()
        
        for filename,data in sorted(self.files.items()): 
            fileItem = QStandardItem((filename))
            fileItem.setDragEnabled(True)
            fileItem.setCheckState(Qt.Unchecked)
            fileItem.setCheckable(True)
            self.series_list_root = self.series_list_model.invisibleRootItem()
            self.series_list_root.appendRow(fileItem)
            self.status_text.setText("Loaded " + filename)
            for name in data.spectra_names():
                item = QStandardItem(name)
                item.setDragEnabled(True)
                item.setCheckState(Qt.Unchecked)
                item.setCheckable(True)
                fileItem.appendRow(item)
    
    def update_ui(self):
        return
    
    def update_file_checks(self, item):
        for file in range(self.series_list_root.rowCount()):
          if id(item) == id(self.series_list_root.child(file)):
            model_index = self.series_list_root.child(file).index()
            checked = self.series_list_model.data(model_index,Qt.CheckStateRole) == Qt.Checked
            for row in range(self.series_list_root.child(file).rowCount()):
              if checked:
                self.series_list_root.child(file).child(row).setCheckState(Qt.Checked)
              else:
                self.series_list_root.child(file).child(row).setCheckState(Qt.Unchecked)
        self.on_show()

    def on_tree_change(self):
        self.update_file_checks()
        self.on_show()

    def on_show(self):
        self.axes.clear()        
        self.axes.set_autoscale_on(self.autoscale_cb.isChecked())
        self.axes.grid(True)
        
        has_series = False
        
        offset = 0
        stacking_param = 1.25
        ticks = array([])
        labels = []

        for file in range(self.series_list_root.rowCount()):
          max_val = 0
          has_series = False; 
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index,
                Qt.CheckStateRole) == (Qt.Checked)
            name = str(self.series_list_model.data(model_index))
            if checked:
                has_series = True
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                if self.BG_cb.isChecked():
                  m = max(spectrum.data)
                else:
                  m = max(spectrum.nobg())
                if m>max_val:
                  max_val=m
           
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index,
                Qt.CheckStateRole) == (Qt.Checked)
            name = str(self.series_list_model.data(model_index))
            
            if checked:
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                scale = 1.0

                if self.NE_cb.isChecked():
                  if self.BG_cb.isChecked():
                    if self.normalize_cb.isChecked():
                      scale = max_val
                    spectrum.plot_full_summary(scale=scale,axes=self.axes, displayParams=self.param_cb.isChecked(),offset=offset)
                  else:
                    if self.normalize_cb.isChecked():
                      scale = max_val
                    spectrum.plot_full_summary_nobg(scale=scale,axes=self.axes, displayParams=self.param_cb.isChecked(),offset=offset)
                if self.dNE_cb.isChecked():
                  spectrum.plot_sg1(scale=scale,points=5,axes=self.axes,offset=offset)

          if self.stacked_cb.isChecked():
            if self.normalize_cb.isChecked():
              if has_series:
                ticks = r_[ticks, offset + r_[0.0:1.01:0.25]]
                self.axes.axis[offset] = self.axes.new_floating_axis(0,offset)
                self.axes.axis[offset].toggle(ticklabels=False)
                new_labels = r_[0.0:1.01:0.25]
                labels += map(str,new_labels)
                offset += stacking_param
            else:
              if has_series:
                ticks = r_[ticks, offset]
                labels += ['$0$']
                self.axes.axis[offset] = self.axes.new_floating_axis(0,offset)
                self.axes.axis[offset].toggle(ticklabels=False)
                inc = (10**floor(log10(max_val)))
                inc = floor(max_val/inc) * inc /2
                for val in r_[inc:max_val*1.01:inc]:
                  ticks = r_[ticks, offset + val]
                  labels += [r'$%.2f \times 10^{%.0f}$'%(val/10**(floor(log10(val))), floor(log10(val)))]
                offset += max_val*stacking_param
        if self.stacked_cb.isChecked():
          self.axes.set_yticks(ticks)
          self.axes.set_yticklabels(labels)
        
        self.canvas.draw()

    def clear_peaks(self):
        print 'clear_peaks'
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                self.files[filename].get_spectrum(row).clear_peaks()
        self.on_show()

    def clear_bg(self):
        print 'clear_bg'
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                print 'bg_cleared ' + filename 
                self.files[filename].get_spectrum(row).clear_bg()
                self.files[filename].get_spectrum(row).clear_abg()
        self.on_show()
    
    def load_fits(self):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFile)
        dialog.setNameFilter('Fit Files (*.fit);;All Files (*.*)')
        if dialog.exec_():
          filename = dialog.selectedFiles()[0]
        f = open(filename, 'r')
        spec = f.read()
        exec(spec)
        spec = spectrum
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                spectrum.set_spec(spec)
        self.on_show()


    def write_fits(self):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly,True)
        if dialog.exec_():
          directory = dialog.selectedFiles()[0]
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                outfilename = filename + '- region' + str(row) + ' - ' + spectrum.name + '.fit'
                print 'write_fits: ' +directory+ '/'+outfilename 
                f = open(directory+'/'+outfilename, 'w')
                f.write('spectrum = ')
                pprint(spectrum.get_spec(), f)
                f.close()

    def adjust_offset(self):
      offset, ok = QInputDialog.getDouble(self, "QInputDialog.getDouble()", "Amount:", 0, -20, 20, 3)
      if ok:
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                spectrum.offset += offset
      self.on_show()

    def plot_optimization_history(self, spectrum):
      for peak in spectrum.peaks:
        if not hasattr(peak, 'optWin'):
          peak.optimization_window = OptimizationWindow(parent=self)

        peak.optimization_window.update()
         
    def optimize_peaks(self):
        for file in range(self.series_list_root.rowCount()):
          for row in range(self.series_list_root.child(file).rowCount()):
            model_index = self.series_list_root.child(file).child(row).index()
            checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
            if checked:
                filename = self.series_list_root.child(file).text()
                spectrum = self.files[filename].get_spectrum(row)
                spectrum.peaks.optimize_fit(spectrum.E(), spectrum.nobg())
                self.plot_optimization_history(spectrum)
        self.on_show()
      
    def write_summary_csv(self):
      for file in range(self.series_list_root.rowCount()):
        for row in range(self.series_list_root.child(file).rowCount()):
          model_index = self.series_list_root.child(file).child(row).index()
          checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
          if checked:
            filename = self.series_list_root.child(file).text()
            spectrum = self.files[filename].get_spectrum(row)
            outfilename = filename + ' - region' + str(row) + ' - ' + spectrum.name + '.csv'
            print 'write_summary: ' + outfilename
            #f = open(outfilename, 'w')
            out = array([])
            out = c_[spectrum.E(), spectrum.data]
            for peak in spectrum.peaks.peak_list:
              out = c_[out,peak(spectrum.E())]
            out = c_[out,spectrum.bg(spectrum.E(),spectrum.data)]
            out = c_[out,spectrum.peaks(spectrum.E())]
            #f.write(out)
            #f.close()
            savetxt(str(outfilename), out, delimiter=",", fmt="%10.5f")

    def plot_summary_csv(self):
      pass

    def average(self):
      i = 0
      outputname = ''
      for file in range(self.series_list_root.rowCount()):
        for row in range(self.series_list_root.child(file).rowCount()):
          model_index = self.series_list_root.child(file).child(row).index()
          checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
          if checked:
            filename = self.series_list_root.child(file).text()
            spectrum = self.files[filename].get_spectrum(row)
            outputname += str(filename+spectrum.name)
            if i==0:
              summer = copy(spectrum.data)
              E = spectrum.E()
            else:
              summer += spectrum.data

            i += 1

      ave = summer/i
      spectrum = Spectrum()

      spectrum.EE = E
      spectrum.data = ave
      spectrum.name = 'Average'

      self.averageWindow = AverageWindow(spectrum, parent=self)
      self.averageWindow.show()
      out = c_[spectrum.EE,spectrum.data]
      np.savetxt(string.replace(outputname,'/','')+'.csv', out, fmt="%12.6G")
      #figure()
 

    def on_about(self):
        msg = __doc__
        QMessageBox.about(self, "About the demo", msg.strip())

    def modify_fit(self, artist):
        self.parameterDialog = ParameterDialog(artist.fit_object, parent=self)
        self.parameterDialog.show()

    def on_pick(self, event):
        if isinstance(event.artist, matplotlib.text.Annotation):
          self.modify_fit(event.artist)
          return True
        try:
            N = len(event.ind)
        except AttributeError:
            return True
        if not N: return True
        x = event.artist.get_data()[0][event.ind[0]]
        y = event.artist.get_data()[1][event.ind[0]]
        if self.normalize_cb.isChecked():
          y *= max(event.artist.spectrum.data)

        self.axes.annotate("(%0.2f,%0.2f)"%(x,y),
                    xy=(x, y), xycoords='data',
                    xytext=(10.5,10.5), textcoords='offset points',
                    arrowprops=dict(arrowstyle="->", connectionstyle="arc3"),
                    bbox=dict(boxstyle="round", fc="w"))
                                                                                                    
        A_scale = 2.0
        if self.sPeakRadio.isChecked():
            peak_spec = { 
                'name': 's - gl',
                'function': 'gl',
                'penalty_function': 'quad_penalty',
                'variables':        ['A',           '\mu',         '\sigma',   'm'],
                'values': r_[ y,            x,             1.0,        0.5],
                'ranges':           [ r_[0, A_scale*y], x+r_[-10, 10], r_[0.1, 5], r_[0.1, 0.9] ]
                }
            event.artist.spectrum.guess_peak_from_spec(peak_spec)
        elif self.pPeakRadio.isChecked():
            peak_spec = { 
                'name': 'p - gl',
                'function': 'spin_split_gl',
                'penalty_function': 'quad_penalty',
                'variables':        ['A', '\mu', '\sigma', 'm', 'R_A', 'R_{Area}', '\Delta E'],
                'values': r_[ y,   x,     1.0,      0.5, 0.5,   0.5,        2],
                'ranges': [ r_[0, A_scale*y],
                            x+r_[-10, 10],
                            r_[0.1, 5],
                            r_[0.1, 0.9],
                            r_[0.42, 0.58],
                            r_[0.49, 0.51],
                            r_[0.1, 20] ]
                }
            event.artist.spectrum.guess_peak_from_spec(peak_spec)
        elif self.dPeakRadio.isChecked():
            peak_spec = { 
                'name': 'd - gl',
                'function': 'spin_split_gl',
                'penalty_function': 'quad_penalty',
                'variables':        ['A', '\mu', '\sigma', 'm', 'R_A', 'R_{Area}', '\Delta E'],
                'values': r_[ y,   x,     1.0,      0.5, 0.5,   0.5,        2],
                'ranges': [ r_[0, A_scale*y],
                            x+r_[-10, 10],
                            r_[0.1, 5],
                            r_[0.1, 0.9],
                            r_[0.586, 0.746],
                            r_[0.656, 0.676],
                            r_[0.1, 20] ]
                }
            event.artist.spectrum.guess_peak_from_spec(peak_spec)
        elif self.fPeakRadio.isChecked():
            peak_spec = { 
                'name': 'f - gl',
                'function': 'spin_split_gl',
                'penalty_function': 'quad_penalty',
                'variables':        ['A', '\mu', '\sigma', 'm', 'R_A', 'R_{Area}', '\Delta E'],
                'values': r_[ y,   x,     1.0,      0.5, 0.5,   0.5,        4.8],
                'ranges': [ r_[0, A_scale*y],
                            x+r_[-10, 10],
                            r_[0.1, 5],
                            r_[0.1, 0.9],
                            r_[0.67, 0.83],
                            r_[0.74, 0.76],
                            r_[4.75, 4.85] ]
                }
            event.artist.spectrum.guess_peak_from_spec(peak_spec)
        elif self.tougaardRadio.isChecked():
            bg_spec = { 'name': 'tougaard',
                        'function': 'K',
                        'penalty_function': 'quad_penalty',
                        'variables':        ['B',          'D'],
                        'values': r_[ 200,         100],
                        'ranges':           [ r_[1, 4000], r_[1, 2000] ]
                        }
            event.artist.spectrum.guess_bg_from_spec(bg_spec)
        elif self.tougaard3Radio.isChecked():
            bg_spec = { 'name': 'tougaard3',
                        'function': 'K3',
                        'penalty_function': 'quad_penalty',
                        'variables':        ['B',          'C',         'D'],
                        'values': r_[ 200,         100,         100], 
                        'ranges':           [ r_[1, 4000], r_[1, 2000], r_[1, 2000] ]  
                        }
            event.artist.spectrum.guess_bg_from_spec(bg_spec)
        elif self.sloped_arctanRadio.isChecked():
            bg_spec = { 'name': 'arctan',
                        'function': 'fitting_arctan',
                        'penalty_function': 'quad_penalty',
                        'variables':  ['A', 'B', 'u'],
                        'values': r_[ y, 1/100.0, x, 0],
                        'ranges': [ r_[ -A_scale*y, A_scale*y],
                                    r_[0, 1000],
                                    x+r_[-100,100],
                                    r_[-1e9, 1e9] ]
                      }
            event.artist.spectrum.guess_abg_from_spec(bg_spec)
        elif False:#self.sloped_sarctanRadio.isChecked():
            bg_spec = { 'name': 'sloped arctan',
                        'function': 'sloped_arctan',
                        'penalty_function': 'quad_penalty',
                        'variables':  ['A', 'B', 'u', 'm', 'b'],
                        'values': r_[ y, 1/100.0, x, 0, 0],
                        'ranges': [ r_[ -A_scale*y, A_scale*y],
                                    r_[0, 1000],
                                    x+r_[-100,100],
                                    r_[-1e9,1e9],
                                    r_[-1e9,1e9] ]
                      }
            event.artist.spectrum.guess_abg_from_spec(bg_spec)

        self.plot_optimization_history(event.artist.spectrum)   
        self.on_show()

        return True

    def create_main_frame(self):
        self.main_frame = QWidget()
        
        plot_frame = QWidget()
        
        self.dpi = 100
        self.fig = Figure((6.0, 4.0), dpi=self.dpi, facecolor='w', edgecolor='k')
        self.canvas = FigureCanvas(self.fig)
        #self.canvas.setParent(self.main_frame)
        self.cidpress = self.fig.canvas.mpl_connect('pick_event', self.on_pick)

        self.axes = Subplot(self.fig,1,1,1)
        self.fig.add_subplot(self.axes)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        
        log_label = QLabel("Spectra:")
        self.series_list_view = QTreeView()
        self.series_list_view.setModel(self.series_list_model)
        self.series_list_view.setDragDropMode(QAbstractItemView.InternalMove)

        #Checkboxes to control plotting display
        self.param_cb = QCheckBox("Show Parameters")
        self.param_cb.setChecked(False)
        self.param_cb.stateChanged.connect(self.on_show)

        self.autoscale_cb = QCheckBox("&Autoscale")
        self.autoscale_cb.setChecked(True)
        self.autoscale_cb.stateChanged.connect(self.on_show)

        self.normalize_cb = QCheckBox("&Normalize")
        self.normalize_cb.setChecked(False)
        self.normalize_cb.stateChanged.connect(self.on_show)

        self.NE_cb = QCheckBox("&NE(E)")
        self.NE_cb.setChecked(True)
        self.NE_cb.stateChanged.connect(self.on_show)

        self.dNE_cb = QCheckBox("&dNE(E)/dE")
        self.dNE_cb.setChecked(False)
        self.dNE_cb.stateChanged.connect(self.on_show)

        self.BG_cb = QCheckBox("&Show Background")
        self.BG_cb.setChecked(True)
        self.BG_cb.stateChanged.connect(self.on_show)

        self.stacked_cb = QCheckBox("&Stacked")
        self.stacked_cb.setChecked(True)
        self.stacked_cb.stateChanged.connect(self.on_show)

        #plotting display controls layout
        cb_hbox1 = QHBoxLayout()
        cb_hbox2 = QHBoxLayout()
        cb_hbox3 = QHBoxLayout()
        cb_vbox = QVBoxLayout()
        cb_hbox1.addWidget(self.param_cb)
        cb_hbox1.addWidget(self.autoscale_cb)
        cb_hbox2.addWidget(self.NE_cb)
        cb_hbox2.addWidget(self.dNE_cb)
        cb_hbox3.addWidget(self.normalize_cb)
        cb_hbox3.addWidget(self.BG_cb)
        cb_hbox3.addWidget(self.stacked_cb)
        cb_vbox.addLayout(cb_hbox1)
        cb_vbox.addLayout(cb_hbox2)
        cb_vbox.addLayout(cb_hbox3)
        
        #select fit types to add when clicking
        self.peakTypeGroup = QGroupBox("Fit Type to Add")
        self.sPeakRadio = QRadioButton("&s")
        self.pPeakRadio = QRadioButton("&p")
        self.dPeakRadio = QRadioButton("&d")
        self.fPeakRadio = QRadioButton("&f")
        self.tougaardRadio = QRadioButton("&Tougaard")
        self.tougaard3Radio = QRadioButton("&Tougaard3")
        self.sloped_arctanRadio = QRadioButton("&arctan")
        self.tougaardRadio.setChecked(True)

        #layout fitting types
        peak_vbox = QVBoxLayout()
        bg_vbox = QVBoxLayout()
        fit_hbox = QHBoxLayout()
        peak_vbox.addWidget(self.sPeakRadio)
        peak_vbox.addWidget(self.pPeakRadio)
        peak_vbox.addWidget(self.dPeakRadio)
        peak_vbox.addWidget(self.fPeakRadio)
        bg_vbox.addWidget(self.tougaardRadio)
        bg_vbox.addWidget(self.tougaard3Radio)
        bg_vbox.addWidget(self.sloped_arctanRadio)
        bg_vbox.addStretch(1)
        fit_hbox.addLayout(peak_vbox)
        fit_hbox.addLayout(bg_vbox)
        self.peakTypeGroup.setLayout(fit_hbox)

        #action buttons
        self.button_clear_peaks = QPushButton("&Clear Peaks")
        self.button_clear_peaks.clicked.connect(self.clear_peaks)
        self.button_clear_peaks.setShortcut("Ctrl+C")

        self.button_optimize_peaks = QPushButton("Optimize &Peaks")
        self.button_optimize_peaks.clicked.connect(self.optimize_peaks)
        self.button_optimize_peaks.setShortcut("Ctrl+P")

        self.button_clear_bg = QPushButton("Clear &BG")
        self.button_clear_bg.clicked.connect(self.clear_bg)
        self.button_clear_bg.setShortcut("Ctrl+B")

        self.button_adjust_offset = QPushButton("&Adjust Offset")
        self.button_adjust_offset.clicked.connect(self.adjust_offset)
        self.button_adjust_offset.setShortcut("Ctrl+A")

        self.button_write_fits = QPushButton("&Write Fits")
        self.button_write_fits.clicked.connect(self.write_fits)
        self.button_write_fits.setShortcut("Ctrl+W")

        self.button_load_fits = QPushButton("&Load Fits")
        self.button_load_fits.clicked.connect(self.load_fits)
        self.button_load_fits.setShortcut("Ctrl+L")

        self.button_write_summary = QPushButton("&Export Summary")
        self.button_write_summary.clicked.connect(self.write_summary_csv)
        self.button_write_summary.setShortcut("Ctrl+E")

        self.button_average = QPushButton("&Average")
        self.button_average.clicked.connect(self.average)
        self.button_average.setShortcut("Ctrl+A")


        #action buttons layout
        mods_box = QGridLayout()
        mods_box.addWidget(self.button_optimize_peaks,0,0)
        mods_box.addWidget(self.button_clear_peaks,0,1)
        mods_box.addWidget(self.button_clear_bg,0,2)
        mods_box.addWidget(self.button_adjust_offset,1,0)
        mods_box.addWidget(self.button_write_fits,1,1)
        mods_box.addWidget(self.button_load_fits,1,2)
        mods_box.addWidget(self.button_write_summary,3,1)
        mods_box.addWidget(self.button_average,3,2)

        left_vbox = QVBoxLayout()
        left_vbox.addWidget(self.canvas)
        left_vbox.addWidget(self.mpl_toolbar)
        left_widget = QWidget()
        left_widget.setLayout(left_vbox)

        right_vbox = QVBoxLayout()
        right_vbox.addWidget(log_label)
        right_vbox.addWidget(self.series_list_view,stretch=1)
        right_vbox.addLayout(cb_vbox)
        right_vbox.addWidget(self.peakTypeGroup)
        right_vbox.addLayout(mods_box)
        right_widget = QWidget()
        right_widget.setLayout(right_vbox)
        
        self.main_frame = QSplitter()
        self.main_frame.addWidget(left_widget)
        self.main_frame.addWidget(right_widget)

        self.setCentralWidget(self.main_frame)
    
    def create_status_bar(self):
        self.status_text = QLabel("Please load a data file")
        self.statusBar().addWidget(self.status_text, 1)

    def create_menu(self):        
        self.file_menu = self.menuBar().addMenu("&File")
        
        load_action = self.create_action("&Open file",
            shortcut="Ctrl+O", slot=self.load_file, tip="Open a file")
        quit_action = self.create_action("&Quit", slot=self.close, 
            shortcut="Ctrl+Q", tip="Close the application")
        
        self.add_actions(self.file_menu, 
            (load_action, None, quit_action))
            
        self.help_menu = self.menuBar().addMenu("&Help")
        about_action = self.create_action("&About", 
            shortcut='F1', slot=self.on_about, 
            tip='About the demo')
        
        self.add_actions(self.help_menu, (about_action,))

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def create_action(  self, text, slot=None, shortcut=None, 
                        icon=None, tip=None, checkable=False, 
                        signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action


class DataHolder(object):
    """ Just a thin wrapper to hold Spectra
    """
    def __init__(self, filename=None, filefilter=None):
        self.load_from_file(filename,filefilter)
    
    def load_from_file(self, filename=None, filefilter=None):
        self.spectra = []
        self.names = []
        
        if filename:
          if filefilter == 'VAMAS (*.vms)':
            experiment = VAMASExperiment(filename)
            for block in experiment.blocks:
                self.names.append(block.sample_identifier + '-' + block.block_identifier)
                spectrum = Spectrum()
                spectrum.EE = block.abscissa()
                spectrum.name = (block.sample_identifier + '-' + block.block_identifier)
                spectrum.data = block.ordinate(0)/(block.number_of_scans_to_compile_this_block * block.signal_collection_time)
                self.spectra.append(spectrum)
          elif filefilter == 'SSRL (*.dat)':
            experiment = np.genfromtxt(filename,skip_header=37)

            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,2]
            spectrum.name = 'I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,3]
            spectrum.name = 'I1'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)
            
            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,2]/experiment[:,1]
            spectrum.name = 'I1/I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,6]
            spectrum.name = 'ICR'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,7]
            spectrum.name = 'FF1'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)
            
            spectrum = Spectrum()
            spectrum.EE = experiment[:,1]
            spectrum.data = experiment[:,6]/experiment[:,1]
            spectrum.name = 'ICR/I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)
          elif filefilter == 'BL7.0.1.1 XAS (*.txt)':
            experiment = np.genfromtxt(filename,
                                       dtype='f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4,f4',
                                       names=['TimeOfDay', 
                                              'Time', 
                                              'MonoEnergy', 
                                              'BeamCurrent', 
                                              'ShutterStatus', 
                                              'Izero', 
                                              'Counter0', 
                                              'Counter1', 
                                              'Counter2', 
                                              'Counter3', 
                                              'Counter4', 
                                              'Counter5', 
                                              'Counter6', 
                                              'Gate7Out', 
                                              'TempA', 
                                              'TempB', 
                                              'TempC', 
                                              'TempD', 
                                              'ColdCathodeGauge', 
                                              'SREnergy'],
                                       skip_header=14)

            spectrum = Spectrum()
            spectrum.EE = experiment['MonoEnergy']
            spectrum.data = experiment['Izero']
            spectrum.name = 'I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment['MonoEnergy']
            spectrum.data = experiment['Counter1']
            spectrum.name = 'TEY'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment['MonoEnergy']
            spectrum.data = experiment['Counter2']
            spectrum.name = 'TFY'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment['MonoEnergy']
            spectrum.data = experiment['Counter1']/experiment['Izero']
            spectrum.name = 'TEY/I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

            spectrum = Spectrum()
            spectrum.EE = experiment['MonoEnergy']
            spectrum.data = experiment['Counter2']/experiment['Izero']
            spectrum.name = 'TFY/I0'
            self.spectra.append(spectrum)
            self.names.append(spectrum.name)

          elif filefilter == 'SUPER (*)':
            pass



    def spectra_names(self):
        """ Names of the data series
        """
        return self.names
    
    def spectra_count(self):
        return len(self.spectra)

    def get_spectrum(self, index):
        return self.spectra[index]


def main():
    matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
    app = QApplication(sys.argv)
    form = Form()
    QApplication.setStyle(QStyleFactory.create('Plastique'))
    QApplication.setPalette(QApplication.style().standardPalette())
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
    print "Fired up"

 #   app = QApplication(sys.argv)
 #   form = Form()
 #   form.show()
 #   QApplication.setStyle(QStyleFactory.create('Plastique'))
 #   QApplication.setPalette(QApplication.style().standardPalette())
    
    
     
