import sys, os  # , csv, string

from spectra_fitting import Spectrum
from VAMAS import *
from data_formats import *

class SpectraFile(object):
    """ Just a thin wrapper to hold Spectra from an individual file
    """
    def __init__(self, filename=None, filefilter=None):
        self.load_from_file(filename,filefilter)

    def load_from_file(self, filename=None, filefilter=None):
        ## fix reload behavoir to maintain existing Spectrum objects
        self.spectra = []
        self.names = []

        if filename:
            if filefilter == 'VAMAS (*.vms)':
                experiment = VAMASExperiment(filename)
                i=1
                for block in experiment.blocks:    #VAMAS file may contain multiple spectra
                    treename = 'B-' + str(i) + '-' + block.comment
                    treename = treename.replace('Not Specified','')    #DFO
                    treename = treename.strip()
                    spectrum = Spectrum()
                    spectrum.EE = block.abscissa()
                    spectrum.name = os.path.basename(filename) + '-' + treename
                    spectrum.data = block.ordinate(0)/(block.number_of_scans_to_compile_this_block * block.signal_collection_time)
                    spectrum.raw_data = spectrum.data.copy()  #keep a copy
                    if block.analyser_mode in ('FAT','FRR') :
                        spectrum.spec_mode = block.analyser_mode
                    spectrum.pass_energy_ratio = block.analyser_pass_energy_of_retard_ratio_or_mass_resolution

                    self.spectra.append(spectrum)
                    self.names.append(treename)
                    i += 1
            elif filefilter == 'SSRL (*.dat)':
                experiment = np.genfromtxt(filename,skip_header=37)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,2]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,3]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'I1'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,2]/experiment[:,1]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'I1/I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,6]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'ICR'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,7]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'FF1'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment[:,1]
                spectrum.data = experiment[:,6]/experiment[:,1]
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'ICR/I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)
            elif filefilter == 'BL7.0.1.1 XAS (*.txt)':
                experiment = read_bl7011_xas(filename)

                spectrum = Spectrum()
                spectrum.EE = experiment['MonoEnergy']
                spectrum.data = experiment['Izero']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment['MonoEnergy']
                spectrum.data = experiment['Counter1']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'TEY'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment['MonoEnergy']
                spectrum.data = experiment['Counter2']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'TFY'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment['MonoEnergy']
                spectrum.data = experiment['Counter1']/experiment['Izero']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'TEY/I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

                spectrum = Spectrum()
                spectrum.EE = experiment['MonoEnergy']
                spectrum.data = experiment['Counter2']/experiment['Izero']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                spectrum.name = 'TFY/I0'
                self.spectra.append(spectrum)
                self.names.append(spectrum.name)

            elif filefilter == 'AugerScan (*.txt)':
                experiment = read_augerscan(filename)
                spectrum = Spectrum()
                spectrum.EE = -experiment['Energy']
                #print spectrum.EE
                spectrum.data = experiment['Counts']
                spectrum.raw_data = spectrum.data.copy()  #keep a copy
                #print spectrum.data
                spectrum.name = 'Region 1'
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
