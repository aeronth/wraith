import os
import re
from pylab import *
from scipy.optimize import leastsq
from fitting_machinery import *
from pprint import *
from scipy import signal
from functions import *

color_list = ['b', 'r', 'g', 'm', 'c', 'y', 'k' ]
color_index = 0

def next_plot_color():
    global color_list
    global color_index
    color = color_list[color_index%len(color_list)]
    color_index += 1
    #assign colors in order for new plots
    #print 'next_color', color
    return color

class Spectrum:
    """A place to keep a spectrum"""

    def __init__(self):
        """Initialize stuff"""
        self.data = ones(1) #raw or pre-processed (pass energy, dead time correct, etc)
        self.raw_data = ones(1)     #raw data read from file,
        self.EE = ones(1)
        self.offset = 0
        self.scale = 1.0
        self.bg = lambda E, data: zeros(E.size)
        self.abg = lambda E: zeros(E.size)
        self.peaks = Peaks(self)
        self.name = "unnamed"
        self.color = ''
        self.spec_mode = '' #FRR or FAT for XPS/Auger
        self.pass_energy_ratio = 1.0        #pass energy eV for FAT, retard ratio FRR

    def pre_process(self,mode):
        if mode:
            #dead time correction, for omicron auger data
            dead_time = 70e-9 / 7.0 #7 channel detector 70 ns dead time
            self.data = self.raw_data / (1 - dead_time * self.raw_data)

            #normalize counts by spec resolution, Hz/eV
            #in FAT mode pass energy is pass_energy_ratio
            #in FRR mode pass energy is KE / pass_energy_ratio
            spec_dispersion = 0.02  #Omicron SCA per-channel resolution/pass energy
            if self.spec_mode == 'FAT':
                self.data = self.data / (spec_dispersion * self.pass_energy_ratio)
            elif self.spec_mode == 'FRR':
                self.data = (self.pass_energy_ratio / spec_dispersion) * self.data / self.EE
        else:
            self.data = self.raw_data.copy()

    def clear_peaks(self):
        self.peaks = Peaks(self)

    def clear_bg(self):
        self.bg = lambda E, data: zeros(E.size)

    def clear_abg(self):
        self.abg = lambda E: zeros(E.size)

    def get_spec(self):
        spec = {}
        spec['offset'] = self.offset
        spec['bg'] = {}
        spec['peaks'] = {}
        try:
            spec['bg'] = self.bg.get_spec()
        except:
            pass

        try:
            spec['peaks'] = self.peaks.get_spec()
        except:
            pass
        return spec

    def set_spec(self, spec):
        self.clear_abg()
        self.bg = Background(self)
        self.peaks = Peaks(self)
        self.offset = spec['offset']
        self.bg.set_spec(spec['bg'])
        self.peaks.set_spec(spec['peaks'])

    def E(self):
        return self.EE - self.offset

    def write_fits(self):
        ffile = open(self.name + '.py', 'w')
        pprint(self.get_spec(), file)
        ffile.close()

    def crop(self, range):
        """Crop to fit within the energy window given in range"""
        upper = self.EE>=range[0]+self.offset
        self.EE = self.EE[upper]
        self.data = self.data[upper]

        lower = self.EE<=range[1]+self.offset
        self.EE = self.EE[lower]
        self.data = self.data[lower]

        #Savitsky Golay smoothing
    def _sg(self,A,points):
        order = 4
        if points<5:
            order = 2
        kernel = savgol(points,points/2,points/2,0,order)#-ker[points]
        pad_start = array([])
        pad_end = array([])
        for i in range(0,(points/2)):
            pad_start = r_[A[0], pad_start]
            pad_end = r_[pad_end, A[-1]]
        sg = r_[pad_start, A, pad_end]
        sg = signal.fftconvolve(sg, kernel, 'valid')
        return sg

    def sg(self,points):
        return self._sg(self.data,points)

        #Savitsky Golay derivative
    def _sg1(self,dE,A,points):
        order = 4
        if points<5:
            order = 2
        kernel = savgol(points, points/2, points/2,1,order)/dE
        pad_start = array([])
        pad_end = array([])
        for i in range(0,(points/2)):
            pad_start = r_[A[0], pad_start]
            pad_end = r_[pad_end, A[-1]]
        sg1 = r_[pad_start, A, pad_end]
        sg1 = convolve(sg1, kernel, 'valid')
        return sg1

    def sg1(self,points):
        return self._sg1(self.E()[0]-self.E()[1], self.data, points)

    def noabg(self):
        """return spectrum with bg subtracted"""
        return self.data - self.abg(self.E())

    def nobg(self):
        """return spectrum with bg subtracted"""
        return self.data - self.abg(self.E()) - self.bg(self.E(), self.data)

    def residuals(self):
        """return residuals from current fit"""
        spectrum = self.nobg() - self.peaks(self.E())
        return spectrum

    def full_fit(self):
        """return the current fit including peaks and bg"""
        return self.peaks(self.E())+self.bg(self.E(),self.data)

    def add_peak(self, name, variables, values, penalties, f):
        """Add a peak with fitting params, using function f"""
        self.peaks.add_peak(Peak(self, name, variables, values, penalties, f))
        self.peaks.optimize_fit(self.E(), self.nobg())

    def guess_abg_from_spec(self, spec):
        values = copy(spec['values'])
        ranges = copy(spec['ranges'])
        self.guess_abg(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

    def guess_abg(self, name, variables, values, ranges, function, penalty_function):
        penalties = []
        for range in ranges:
            penalties.append(Penalty(range,penalty_function))

        self.abg = AnalyticBackground(self, name, variables, values, penalties, function)
        self.abg.optimize_fit(self.E(), self.data)

    def guess_bg_from_spec(self, spec):
        values = copy(spec['values'])
        ranges = copy(spec['ranges'])
        self.guess_bg(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

    def guess_bg(self, name, variables, values, ranges, function, penalty_function):
        penalties = []
        for range in ranges:
            penalties.append(Penalty(range,penalty_function))

        self.bg = Background(self, name, variables, values, penalties, function, values[1]*2)
        self.bg.optimize_fit(self.E(), self.noabg())

    def guess_peak_from_spec(self, spec):
        values = copy(spec['values'])
        ranges = copy(spec['ranges'])

        self.guess_peak(spec['name'], spec['variables'], values, ranges, eval(spec['function']), eval(spec['penalty_function']))

    def guess_peak(self, name, variables, values, ranges, function, penalty_function):
        penalties = []
        for range in ranges:
            penalties.append(Penalty(range,penalty_function))

        self.add_peak(name, variables, values, penalties, function)

    #various curve and fit plotting routines
    def plot_f_of_data(self,f,scale=1.0,axes=None, offset=0.0):
        #this seems to be left over not called - DFO
        if axes==None:
            axes=gca()
        lines = axes.plot(self.E(), offset+f((self.data)/scale),'--',label=f.func_name + "(" + self.name + ")")

    def plot(self,scale=1.0,axes=None, offset=0.0):
        """Plot the spectra"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset+(self.data)/scale, markersize=3, color=self.color, label=self.name,picker=2)
        for line in lines:
            line.spectrum = self


    def plot_sg(self,points=5, scale=1.0,axes=None, offset=0.0):
        """Plot the smoothed spectra"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset+(self.sg(points))/scale, markersize=3, color=self.color, label=self.name,picker=2)
        for line in lines:
            line.spectrum = self

    def plot_sg1(self,points=5,scale=1.0,axes=None, offset=0.0):
        """Plot the Savitski Golay derviative of the spectra"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset+(self.sg1(points))/scale, color = self.color, markersize=3, label=self.name+'-sg1(%d)'%(points,))

    def plot_nobg(self,scale=1.0, axes=None, offset=0.0):
        """Plot spectrum with bg subtracted"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset+(self.nobg())/scale, markersize=3, color=self.color, label=self.name+'-nobg',picker=2)
        for line in lines:
            line.spectrum = self

    def plot_sg_nobg(self,points=5, scale=1.0,axes=None, offset=0.0):
        """Plot the smoothed spectra"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset+(self._sg(self.nobg(),points))/scale, color=self.color, markersize=3, label=self.name,picker=2)
        for line in lines:
            line.spectrum = self

    def plot_peaks(self,scale=1.0, axes=None, offset=0.0):
        """Plot the peaks summed together"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset + (self.peaks(self.E()))/scale, ':', label=self.name+'-peaks')
        for line in lines:
            line.spectrum = self

    def plot_individual_peaks(self,scale=1.0, axes=None, offset=0.0):
        """Plot the peaks as individual functions"""
        if axes==None:
            axes = gca()
        i = 0
        colors = ['b','g','r','c','m','y','k']
        for peak in self.peaks.peak_list:
            line = axes.fill_between(self.E(),offset, offset + (peak(self.E()))/scale,label=self.name+'-peak%d'%i,alpha=0.4,color=colors[i%len(colors)])
            line.spectrum = self
            peak.line = line
            i += 1

    def plot_individual_peaks_bg(self,scale=1.0, axes=None, offset=0.0):
        """Plot the peaks as individual functions"""
        if axes==None:
            axes = gca()
        b = self.bg(self.E(), self.data)
        i = 0
        colors = ['b','g','r','c','m','y','k']
        for peak in self.peaks.peak_list:
            line = axes.fill_between(self.E(), b/scale+offset, offset + (peak(self.E())+b)/scale,label=self.name+'-peak%d'%i,alpha=0.4,color=colors[i])
            line.spectrum = self
            peak.line = line
            i += 1

    def plot_residuals(self,scale=1.0, axes=None, offset=0.0):
        """plot the residuals"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset + (self.residuals())/scale,label=self.name+'-residuals')
        for line in lines:
            line.spectrum = self

    def plot_bg(self,scale=1.0, axes=None, offset=0.0):
        """Plot the bg"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(), offset + (self.abg(self.E()) + self.bg(self.E(), self.data))/scale,label=self.name+'-bg')
        for line in lines:
            line.spectrum = self
        self.bg.line = lines[0]

    def plot_abg(self,scale=1.0, axes=None, offset=0.0):
        """Plot the bg"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(),offset+(self.abg(self.E()))/scale,label=self.name+'-abg')
        for line in lines:
            line.spectrum = self

    def plot_kbg(self,scale=1.0, axes=None, offset=0.0):
        """Plot the bg"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(),offset+(self.bg(self.E(), self.data))/scale,label=self.name+'-kbg')
        for line in lines:
            line.spectrum = self

    def plot_full_fit(self,scale=1.0, axes=None, offset=0.0):
        """Plot the current fit including peaks and bg"""
        if axes==None:
            axes = gca()
        lines = axes.plot(self.E(),offset + (self.peaks(self.E())+self.bg(self.E(),self.data))/scale,label=self.name+'-fullfit')
        for line in lines:
            line.spectrum = self

    def plot_full_summary(self,smoothpoints=1,scale=1.0, axes=None, displayParams=True, offset=0.0):
        """Plot the spectrum, the bg, and the full fit"""
        if axes==None:
            axes = gca()
        #if smoothpoints<5:
        #  self.plot(scale, axes, offset)
        #else:
        self.plot_sg(smoothpoints, scale, axes, offset)
        self.plot_bg(scale, axes, offset)
        self.plot_abg(scale, axes, offset)
        self.plot_kbg(scale, axes, offset)
        self.plot_full_fit(scale, axes, offset)
        self.plot_individual_peaks_bg(scale, axes, offset)

        if displayParams:
            if isinstance(self.bg, Background):
                spec = self.bg.get_spec()
                label = spec['name'] + "\n"
                label += spec['function']
                for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
                    label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
                label += ""
                mu = self.E()[0]
                A = self.bg(self.E(),self.data)[0]
                EE = self.E()[::-1]
                col = self.bg.line.get_color()
                an = axes.annotate(label,
                              xy=(mu, A), xycoords='data',
                              xytext=(0,60), textcoords='offset points',
                              arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0",fc='%s'%col,alpha=0.6),
                              bbox=dict(boxstyle="roundtooth",fc='%s'%col,alpha=0.2),
                              fontsize='small',
                              family='sans-serif')
                an.set_picker(True)
                #an.draggable()
                an.fit_object = self.bg

            for peak in self.peaks:
                spec = peak.get_spec()
                label = spec['name'] + '\n'
                label += spec['function']
                for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
                    label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
                label += ""
                mu = spec['values'][1]
                A = spec['values'][0]/scale
                EE = self.E()[::-1]
                back = self.bg(self.E(), self.data)[::-1]
                if mu > EE[-1]:
                    b = back[-1]
                else:
                    b = back[EE>=mu][0]
                col = peak.line.get_facecolor()
                an = axes.annotate(label,
                              xy=(mu, A + b), xycoords='data',
                              xytext=(0,60), textcoords='offset points',
                              arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0"),
                              bbox=dict(boxstyle="round"),
                              fontsize='small',
                              family='sans-serif')
                an.get_bbox_patch().set_facecolor(col[0])
                an.get_bbox_patch().set_alpha(col[0][3])
                an.set_picker(True)
                #an.draggable()
                an.fit_object = peak

    def plot_full_summary_nobg(self,smoothpoints=1,scale=1.0, axes=None, displayParams=True, offset=0.0):
        """Plot the spectrum, the bg, and the full fit"""
        if axes==None:
            axes = gca()
        if smoothpoints<5:
            self.plot_nobg(scale, axes, offset)
        else:
            self.plot_sg_nobg(smoothpoints, scale, axes, offset)
        self.plot_peaks(scale, axes, offset)
        self.plot_individual_peaks(scale, axes, offset)

        if displayParams:
            for peak in self.peaks:
                spec = peak.get_spec()
                label = spec['name'] + '\n'
                label += spec['function']
                for var, val, range in zip(spec['variables'], spec['values'], spec['ranges']):
                    label += "\n$%s = %0.2f\in\, (%0.2f,%0.2f)$"%(var,val,range[0],range[1])
                label += ""
                mu = spec['values'][1]
                A = spec['values'][0]/scale
                EE = self.E()[::-1]
                col = peak.line.get_facecolor()
                an = axes.annotate(label,
                              xy=(mu, A), xycoords='data',
                              xytext=(0,60), textcoords='offset points',
                              arrowprops=dict(arrowstyle='fancy',connectionstyle="arc3,rad=0.0"),
                              bbox=dict(boxstyle="round"),
                              fontsize='small',
                              family='sans-serif')
                an.get_bbox_patch().set_facecolor(col[0])
                an.get_bbox_patch().set_alpha(col[0][3])
                an.set_picker(True)
                #an.draggable()
                an.fit_object = peak

    def __call__(self):
        return self.data
