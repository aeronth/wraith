#Main window to control plotting
from PySide.QtCore import *
from PySide.QtGui import *
import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure

from pylab import *
import os

from spectrum_handling import SpectraFile
from data_formats import *
import mpl_toolkits.axisartist as AA
from optimization_gui import *
from parameter_gui import *
from spectra_fitting import *


class Form(QMainWindow):
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)
        self.setWindowTitle('Wraith Interactive Spectral Explorer')
        self.ignore_signals = False

        self.files = {}    #list of open file basenames
        self.series_list_model = QStandardItemModel()
        self.series_list_root = self.series_list_model.invisibleRootItem()
        self.series_list_model.itemChanged.connect(self.update_file_checks)

        self.create_menu()
        self.create_main_frame()
        self.create_status_bar()
        #self.update_ui()
        self.on_show()

    def filterSelected(self, typefilter):    #file type filter set in UI file dialog
        self.filefilter = typefilter

    def clear_files(self):  #clear all files/plots without restarting app, FIX not quite right
        self.files = {}
        self.series_list_model.clear()
        self.series_list_root = self.series_list_model.invisibleRootItem()  #update after clear??
        #self.on_show()  #clear tree control
        self.canvas.draw()  #clear plot area
    
    def load_file(self, filename=None):
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.ExistingFiles)
        dialog.setNameFilter('VAMAS (*.vms);; AugerScan (*.txt);; BL7.0.1.1 XAS (*.txt);; SSRL (*.dat);; SUPER (*);; All Files (*.*)')
        self.filefilter = 'VAMAS (*.vms)'
        dialog.filterSelected.connect(self.filterSelected)

        if dialog.exec_():  #FIX failure to add additional file should not kill list
            filelist = dialog.selectedFiles()
        else:
            filelist = []

        for filename in filelist:
            if filename:
                name = os.path.basename(str(filename))
                self.files[name] = SpectraFile()
                self.files[name].load_from_file(str(filename), str(self.filefilter))
                self.fill_series_list()
                #self.update_ui()

    def fill_series_list(self):
        self.series_list_model.clear()

        for filename,data in sorted(self.files.items()):
            fileItem = QStandardItem((filename))
            fileItem.setDragEnabled(False)    #no functionality yet
            fileItem.setCheckState(Qt.Unchecked)
            fileItem.setCheckable(True)
            self.series_list_root = self.series_list_model.invisibleRootItem()
            self.series_list_root.appendRow(fileItem)
            self.status_text.setText("Loaded " + filename)
            for name in data.spectra_names():
                item = QStandardItem(name)
                item.setDragEnabled(False)
                item.setCheckState(Qt.Unchecked)
                item.setCheckable(True)
                fileItem.appendRow(item)

    #def update_ui(self):
    #    return

    def update_file_checks(self, item):
        for file in range(self.series_list_root.rowCount()):
            if id(item) == id(self.series_list_root.child(file)):
                model_index = self.series_list_root.child(file).index()
                checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == Qt.Checked
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
        if self.ignore_signals:
            return

        self.axes.clear()
        self.axes.set_autoscale_on(self.autoscale_cb.isChecked())
        self.axes.grid(True)

        has_series = False

        offset = 0
        stacking_param = 1.25
        ticks = array([])
        labels = []

        sg_points = 2 * self.smooth_points.value() + 1

        #logic to normalize data by file, may want to do by spectrum instead...DFO
        #the first time through data is not yet loaded..??..DFO
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
                        if self.NE_cb.isChecked():
                            m = max(spectrum.data)
                        else:
                            m = max(spectrum.sg1(sg_points))
                    else:
                        m = max(spectrum.nobg())
                    if m>max_val:
                        max_val=m
            #loop over lines and do plots
            for row in range(self.series_list_root.child(file).rowCount()):
                model_index = self.series_list_root.child(file).child(row).index()
                checked = self.series_list_model.data(model_index,
                    Qt.CheckStateRole) == (Qt.Checked)
                name = str(self.series_list_model.data(model_index))

                if checked:
                    filename = self.series_list_root.child(file).text()
                    spectrum = self.files[filename].get_spectrum(row)
                    spectrum.pre_process(self.pass_energy_cb.isChecked())

                    scale = 1.0
                        #assign color first time, keep for session
                    if spectrum.color == '':
                        spectrum.color = next_plot_color()

                    if self.NE_cb.isChecked():
                        if self.BG_cb.isChecked():
                            if self.normalize_cb.isChecked():
                                scale = max_val
                            spectrum.plot_full_summary(scale=scale, smoothpoints=sg_points, axes=self.axes, displayParams=self.param_cb.isChecked(),offset=offset)
                        else:
                            if self.normalize_cb.isChecked():
                                scale = max_val
                            spectrum.plot_full_summary_nobg(scale=scale, smoothpoints=sg_points, axes=self.axes, displayParams=self.param_cb.isChecked(),offset=offset)
                    if self.dNE_cb.isChecked():

                        spectrum.plot_sg1(scale=scale, points=sg_points, axes=self.axes, offset=offset)
                        self.status_text.setText(spectrum.name)

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
        for specfile in range(self.series_list_root.rowCount()):
            for row in range(self.series_list_root.child(specfile).rowCount()):
                model_index = self.series_list_root.child(specfile).child(row).index()
                checked = self.series_list_model.data(model_index, Qt.CheckStateRole) == (Qt.Checked)
                if checked:
                    filename = self.series_list_root.child(specfile).text()
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
                    #print 'bg_cleared ' + filename
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
#        spec = spectrum
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
        print "plot opt history"
        if not hasattr(spectrum.bg, 'optimization_window'):
            spectrum.bg.optimization_window = OptimizationWindow(spectrum.bg, parent=self)
        spectrum.bg.optimization_window.update()

        for peak in spectrum.peaks:
            if not hasattr(peak, 'optimization_window'):
                peak.optimization_window = OptimizationWindow(peak, parent=self)

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
                    #FIX error here if no peaks to fit...
                    #self.plot_optimization_history(spectrum)
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
            #change to physical def DFO
            bg_spec = { 'name': 'Tougaard',
                        'function': 'K',
                        'penalty_function': 'quad_penalty',
                        'variables':        ['R_loss','E_loss'],
                        'values': r_[ .872,         63.6],
                        'ranges':           [ r_[0.8, 1.0], r_[10.0, 70.0] ]
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

        #self.plot_optimization_history(event.artist.spectrum)
        self.on_show()

        return True

    def create_main_frame(self):
        self.main_frame = QWidget()

        plot_frame = QWidget()

        #set up matplotlib plotting window
        self.dpi = 100
        self.fig = Figure((6.0, 4.0), dpi=self.dpi, facecolor='w', edgecolor='k')
        self.canvas = FigureCanvas(self.fig)
        #self.canvas.setParent(self.main_frame)
        self.cidpress = self.fig.canvas.mpl_connect('pick_event', self.on_pick)

        self.axes = AA.Subplot(self.fig,1,1,1)
        self.fig.add_subplot(self.axes)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)

        #create tree view for data files/spectra
        log_label = QLabel("Spectra:")
        self.treeView = QTreeView()
        self.treeView.setModel(self.series_list_model)
        self.treeView.setDragDropMode(QAbstractItemView.InternalMove)
        self.treeView.setUniformRowHeights(True)
        self.treeView.setHeaderHidden(True)

        #Checkboxes to control plotting display
        #add pass energy correction for Auger, change defaults DFO
        self.param_cb = QCheckBox("Show Param")
        self.param_cb.setChecked(False)
        self.param_cb.stateChanged.connect(self.on_show)

        self.autoscale_cb = QCheckBox("&Autoscale")
        self.autoscale_cb.setChecked(True)
        self.autoscale_cb.stateChanged.connect(self.on_show)

        self.pass_energy_cb = QCheckBox("&Pass Energy")
        self.pass_energy_cb.setChecked(True)
        self.pass_energy_cb.stateChanged.connect(self.on_show)

        self.normalize_cb = QCheckBox("&Normalize")
        self.normalize_cb.setChecked(False)
        self.normalize_cb.stateChanged.connect(self.on_show)

        self.NE_cb = QCheckBox("&NE(E)")
        self.NE_cb.setChecked(True)
        self.NE_cb.stateChanged.connect(self.on_show)

        self.dNE_cb = QCheckBox("&dNE(E)/dE")
        self.dNE_cb.setChecked(False)
        self.dNE_cb.stateChanged.connect(self.on_show)

        self.smooth_points = QSpinBox()    #SG smooth/diff
        self.smooth_points.setRange(1,35)
        self.smooth_points.setValue(1)
        self.smooth_points.valueChanged.connect(self.on_show)

        self.BG_cb = QCheckBox("&Show Background")
        self.BG_cb.setChecked(False)
        self.BG_cb.stateChanged.connect(self.on_show)

        self.stacked_cb = QCheckBox("&Stacked")
        self.stacked_cb.setChecked(False)
        self.stacked_cb.stateChanged.connect(self.on_show)

        #plotting display controls layout
        cb_hbox1 = QHBoxLayout()
        cb_hbox2 = QHBoxLayout()
        cb_hbox3 = QHBoxLayout()
        cb_vbox = QVBoxLayout()
        cb_hbox1.addWidget(self.param_cb)
        cb_hbox1.addWidget(self.autoscale_cb)
        cb_hbox1.addWidget(self.pass_energy_cb)
        cb_hbox2.addWidget(self.NE_cb)
        cb_hbox2.addWidget(self.dNE_cb)
        cb_hbox2.addWidget(self.smooth_points)
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

        #action buttons layout
        mods_box = QGridLayout()
        mods_box.addWidget(self.button_optimize_peaks,0,0)
        mods_box.addWidget(self.button_clear_peaks,0,1)
        mods_box.addWidget(self.button_clear_bg,0,2)
        mods_box.addWidget(self.button_adjust_offset,1,0)
        mods_box.addWidget(self.button_write_fits,1,1)
        mods_box.addWidget(self.button_load_fits,1,2)
        mods_box.addWidget(self.button_write_summary,3,1)

        left_vbox = QVBoxLayout()
        left_vbox.addWidget(self.canvas)
        left_vbox.addWidget(self.mpl_toolbar)
        left_widget = QWidget()
        left_widget.setLayout(left_vbox)

        right_vbox = QVBoxLayout()
        right_vbox.addWidget(log_label)
        right_vbox.addWidget(self.treeView,stretch=1)
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
        load_action = self.create_action("&Open file", slot=self.load_file,
            shortcut="Ctrl+O", tip="Open a file")
        clear_action = self.create_action("&Clear files", slot=self.clear_files,
            shortcut="Ctrl+C", tip="Remove all files and associated plots")
        quit_action = self.create_action("&Quit", slot=self.close,
            shortcut="Ctrl+Q", tip="Close the application")

        self.add_actions(self.file_menu, (load_action, clear_action, quit_action))

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
            #note ToolTips do not work for Qt menu actions (they do for menus)
            #the StatusTips DO work
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            self.connect(action, SIGNAL(signal), slot)
        if checkable:
            action.setCheckable(True)
        return action
