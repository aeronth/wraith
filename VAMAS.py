
#######################################################################################################
# Copyright (c) 2010
#
# Aaron Tynes Hammack
#   Molecular Foundary
#   Lawrence Berkeley National Laboratories
#
# Python module to read in VAMAS files for use with SciPy.
#
# Data import and export proceeds according to the specifications in
#
#   'VAMAS Surface Chemical Analysis Standard Data Transfer Format with Skeleton Decoding Programs'
#       W.A. Dench, L.B. Hazell, M.P. Seah, and the VAMAS Community
#       Surface and Interface Analysis, Volume 13, pages 63-122 (1988)
#
#######################################################################################################

# A note on units standardized in the VAMAS 1988 specification:
# units = ( 'c/s' | 'd' | 'degree' | 'eV' | 'K' | 'micro C' | 'micro m' | 'm/s' | 'n' | 'nA' | 'ps' | 's' | 'u' | 'V')
# These values are abbreviations for the units listed below:
#  'c/s'      counts per second
#  'd'        dimensionless - just a number, e.g. counts per channel
#  'degree'   angle in degrees
#  'eV'       eletron volts
#  'K'        Kelvin
#  'micro C'  microcoulombs
#  'micro m'  micrometres
#  'm/s'      metres per second
#  'n'        not defined here - may be given in a lobel
#  'nA'       nanoamps
#  'ps'       picoseconds
#  's'        seconds
#  'u'        unified atomic mass units
#  'V'        volts


# Import various modules for data and processing
from datetime import datetime
from scipy import *

#Class to store the experiment data contained in a VAMAS file
class VAMASExperiment:
    def set_defaults(self):
        self.format_identifier = "VAMAS Surface Chemical Analysis Standard Data Transfer Format 1988 May 4"
        self.institution_identitifier = "Not Specified"
        self.instrument_model_identifier = "Not Specified"
        self.operator_identifier = "Not Specified"
        self.experiment_identifier = "Not Specified"
        self.number_of_lines_in_comment = 0
        self.comment = "Not Specified"

        self.experiment_mode = "NORM"
        # ( 'MAP' | 'MAPDP' | 'MAPSV' | 'MAPSVDP' | 'NORM' | 'SDP' | 'SDPSV' | 'SEM' ),
        # The contents of each block in the experiment are indicated by the values of experiment mode as follows:
        # 'MAP'
        #     A spectrum which refers to a specified point in a regular two-dimensional spatial array.
        # 'MAPDP'
        #     A spectrum which refers to a specified point in a regular two-dimensional spatial array
        #       and to a specified layer in a depth profile.
        # 'MAPSV'
        #     A complete set of single values of a fixed number of variables for every point in a
        #     regular two-dimensional spatial array. Note that an x linescan consists of a map with
        #     the value of number_of_analysis_positions equal to the value of
        #     number_of_discrete_x_coordinates_available_in_full_map, that is, the number of discrete
        #     y coordinates is unity; in a y linescan the roles of x and y are reversed.
        # 'MAPSVDP'
        #     A complete set of single values of a fixed number of variables for every point in a
        #     regular two-dimensional array for one layer in a depth profile. Successive blocks
        #     refer to successive layers in the depth profile.
        # 'NORM'
        #     Either independent data or data which refers to a specified set of single values of
        #     one or more experimental variables; the data may be spectral or non-spectral.
        # 'SDP'
        #     A spectrum which refers to a specified layer in a depth profile.
        # 'SDPSV'
        #     A complete set of single values of a fixed number of variables for every layer in a depth profile.
        # 'SEM'
        #     An electron emission intensity for every point in a regular two-dimensional spatial array.

        self.scan_mode = "REGULAR"
        # ( 'REGULAR' | 'IRREGULAR' | 'MAPPING' ),
        # If the value of experiment mode is 'MAPSV', 'MAPSVDP' or 'SEM' then the value of scan_mode
        # must be 'MAPPING', otherwise if the data is in the form of an abscissa start, an abscissa
        # increment and a number of complete sets of values of one or more experimental variables then
        # the value of scan mode is 'REGULAR', otherwise the value of scan mode is 'IRREGULAR'.

        # only present if experiement_mode is 'MAP', 'MAPDP', 'NORM', or 'SDP'
        self.number_of_spectral_regions = 0

        # only present if experiement_mode is 'MAP', or 'MAPDP'
        self.number_of_analysis_positions = 0
        self.number_of_discrete_x_coordinates_in_full_map = 0
        self.number_of_discrete_y_coordinates_in_full_map = 0

        # An experimental variable is a parameter which may be varied from block to block
        # through the experiment but which remains constant within each block
        self.number_of_experimental_variables = 0

        # The number of occurances of experiment_variable_labels and _units is specified by
        # the value of number_of_experiment_variables
        self.experimental_variable_labels = []
        self.experimental_variable_units = []

        # If the values of any of the block parameters are the same in all of the blocks
        # their values may be sent in the first block and then omitted from all subsequent
        # blocks. If number_of_entries_in_parameter_inclusion_or_exclusion_list is positive
        # then the parameters listed are to be included; if it is negative then the
        # parameters listed are to be excluded; if it is zero then all parameters are to
        # be given in all blocks.
        self.number_of_entries_in_parameter_inclusion_or_exclusion_list = 0

        # The number of occurrences of parameter_inclusion_or_exclusion_prefix_numbers is given by the
        # absolute value of number_of_entries_in_parameter_inclusion_or_exclusion_list above. If this
        # is greater than zero then the values of successive occurrences of
        # parameter_inclusion_or_exclusion_prefix_numbers should be in ascending order. The values of
        # parameter_inclusion_or_exclusion_prefix_numbers refer to the numbers in comment brackets given below to the left
        # of the parameters in the syntax-rule defining block and must be within the range of these numbers.
        self.parameter_inclusion_or_exclusion_prefix_numbers = []

        self.number_of_manually_entered_items_in_block = 0

        # The number of occurrences of prefix_number_of_manually_entered_items is specified by the value of
        # number of manually entered items in block above. If this is greater than zero then the values of
        # successive occurrences of prefuc number of manually entered item should be in ascending order.
        # Any of the items preceded by prefix numbers in comment brackets in the syntax-rule defining block
        # which need to be evaluated by the operator and manually entered from the keyboard should be included
        # in this list. If an item is to be expressed as a real number and the operator is unable to supply
        # a value then the computer should enter the value 1E37.
        self.prefix_numbers_of_manually_entered_items = []

        # number of future upgrade experiment entries and number of future upgrade block entries are included
        # in case the Format is upgraded in the future to include more non-optional, non-repeating parameters.
        # The numbers of these new parameters will be entered here so that old programs can skip the new
        # parameters in new data, and new programs will not try to read the new parameters in old data.For the
        # present both of them would be set to zero.
        self.number_of_future_upgrade_experiment_entries = 0
        self.number_of_future_upgrade_block_entries = 0

        # The number of occurrences of future upgrade experiment entry is given by the value of number of
        # future upgrade experiment entries above. It is defined as a text line so that any integer,
        # real number or text line inserted here by a future upgrade of the Format can be read as a text line
        # then discarded.
        self.future_upgrade_experiment_entries = []

        self.number_of_blocks = 1

        self.blocks = [] #Array of the VAMASBlocks in the file

        self.experiment_terminator = 'end of experiment'
    #end def __init__(self):

    def __init__(self,filename):
        self.set_defaults()
        self.file = open(filename)
        self.read_from_file()

    def read_from_file(self):
        self.file.seek(0)
        lines = iter(self.file.readlines())

        self.format_identifier = lines.next().strip()
        self.institution_identifier = lines.next().strip()
        self.instrument_model_identifier = lines.next().strip()
        self.operator_identifier = lines.next().strip()
        self.experiment_identifier = lines.next().strip()
        self.number_of_lines_in_comment = int(lines.next().strip())

        self.comment = ""
        for i in range(self.number_of_lines_in_comment):
            self.comment = self.comment + lines.next()

        self.experiment_mode = lines.next().strip()
        self.scan_mode = lines.next().strip()

        if (self.experiment_mode.upper() == 'MAP' or
            self.experiment_mode.upper() == 'MAPDP' or
            self.experiment_mode.upper() == 'NORM' or
            self.experiment_mode.upper() == 'SDP'):
            self.number_of_spectral_regions = int(lines.next().strip())

        if (self.experiment_mode.upper() == 'MAP' or
            self.experiment_mode.upper() == 'MAPDP'):
            self.number_of_analysis_positions = int(lines.next().strip())
            self.number_of_discrete_x_coordinates_available_in_full_map = int(lines.next().strip())
            self.number_of_discrete_y_coordinates_available_in_full_map = int(lines.next().strip())

        self.number_of_experimental_variables = int(lines.next().strip())

        for i in range(self.number_of_experimental_variables):
            self.experimental_variable_labels.append(lines.next().strip())
            self.experimental_variable_units.append(lines.next().strip())

        self.number_of_entries_in_parameter_inclusion_or_exclusion_list = int(lines.next().strip())

        for i in range(abs(self.number_of_entries_in_parameter_inclusion_or_exclusion_list)):
            self.parameter_inclusion_or_exclusion_prefix_numbers.append(int(lines.next().strip()))

        self.number_of_manually_entered_items_in_block = int(lines.next().strip())

        for i in range(self.number_of_manually_entered_items_in_block):
            self.prefix_numbers_of_manually_entered_items.append(int(lines.next().strip()))

        self.number_of_future_upgrade_experiment_entries = int(lines.next().strip())

        self.number_of_future_upgrade_block_entries = int(lines.next().strip())

        for i in range(self.number_of_future_upgrade_block_entries):
            self.future_upgrade_experiment_entries.append(lines.next().strip())

        self.number_of_blocks = int(lines.next().strip())

        for i in range(self.number_of_blocks):
            block = VAMASBlock()
            block.block_identifier = lines.next().strip()
            block.sample_identifier = lines.next().strip()
            # 1
            year = int(lines.next().strip())

            # 2
            month = int(lines.next().strip())

            # 3
            day = int(lines.next().strip())

            # 4
            hours = int(lines.next().strip())
            if hours == 24:
                hours = 0 #catch a 24:00 hours ambiguity

            # 5
            minutes = int(lines.next().strip())

            # 6
            seconds = int(lines.next().strip())
            try:
                block.date = datetime(year,month,day,hours,minutes,seconds)
            except:
                print "odd date: %dy, %dm, %dd, %dh, %dm, %ds"%(year,month,day,hours,minutes,seconds)
                block.date = datetime(1970,1,1,1,1,1)
            # 7
            block.number_of_hours_in_advance_of_greenwich_mean_time = int(lines.next().strip())

            # 8
            block.number_of_lines_in_block_comment = int(lines.next().strip())
            for i in range(block.number_of_lines_in_block_comment):
                block.comment = block.comment + lines.next()

            # 9
            block.technique = lines.next().strip()

            # 10
            if (self.experiment_mode.upper() == 'MAP' or
                self.experiment_mode.upper() == 'MAPDP'):
                block.x_coordinate = int(lines.next().strip())
                block.y_coordinate = int(lines.next().strip())

            # 11
            for i in range(self.number_of_experimental_variables):
                block.values_of_experimental_variables.append(float(lines.next().strip()))

            # 12
            block.analysis_source_label = lines.next().strip()

            # 13
            if (self.experiment_mode.upper() == 'MAPDP' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SDP' or
                self.experiment_mode.upper() == 'SVDP' or
                block.technique.upper() == 'FABMS' or
                block.technique.upper() == 'FABMS energy spec'.upper() or
                block.technique.upper() == 'ISS' or
                block.technique.upper() == 'SIMS' or
                block.technique.upper() == 'SIMS energy spec'.upper() or
                block.technique.upper() == 'SNMS' or
                block.technique.upper() == 'SNMS energy spec'.upper()):
                block.sputtering_ion_or_atomic_number = lines.next().strip()
                block.number_of_atoms_in_sputtering_ion_or_atom_particle = int(lines.next().strip())
                block.sputtering_ion_or_atom_charge_sign_and_number = int(lines.next().strip())

            # 14
            block.analysis_source_characteristic_energy = float(lines.next().strip())

            # 15
            block.analysis_source_strength = float(lines.next().strip())

            # 16
            block.analysis_source_beam_width_x = float(lines.next().strip())
            block.analysis_source_beam_width_y = float(lines.next().strip())

            # 17
            if (self.experiment_mode.upper() == 'MAP' or
                self.experiment_mode.upper() == 'MAPDP' or
                self.experiment_mode.upper() == 'MAPSV' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SEM'):
                block.field_of_view_x = float(lines.next().strip())
                block.field_of_view_y = float(lines.next().strip())

            # 18
            if (self.experiment_mode.upper() == 'MAPSV' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SEM'):
                block.first_linescan_start_x_coordinate = int(lines.next().strip())
                block.first_linescan_start_y_coordinate = int(lines.next().strip())
                block.first_linescan_finish_x_coordinate = int(lines.next().strip())
                block.first_linescan_finish_y_coordinate = int(lines.next().strip())
                block.last_linescan_finish_x_coordinate = int(lines.next().strip())
                block.last_linescan_finish_y_coordinate = int(lines.next().strip())

            # 19
            block.analysis_source_polar_angle_of_incidence = float(lines.next().strip())

            # 20
            block.analysis_source_azimuth = float(lines.next().strip())

            # 21
            block.analyser_mode = lines.next().strip()

            # 22
            block.analyser_pass_energy_of_retard_ratio_or_mass_resolution = float(lines.next().strip())

            # 23
            if (block.technique.upper() == 'AES DIFF'):
                block.differential_width = float(lines.next().strip())

            # 24
            block.magnification_of_analyser_transfer_lens = float(lines.next().strip())

            # 25
            block.analyser_work_function_or_acceptance_energy_of_atom_or_ion = float(lines.next().strip())

            # 26
            block.target_bias = float(lines.next().strip())

            # 27
            block.analysis_width_x = float(lines.next().strip())
            block.analysis_width_y = float(lines.next().strip())

            # 28
            block.analyser_axis_take_off_polar_angle = float(lines.next().strip())
            block.analyser_axis_take_off_azimuth = float(lines.next().strip())

            # 29
            block.species_label = lines.next().strip()

            # 30
            block.transition_or_charge_state_label = lines.next().strip()
            block.charge_of_detected_particle = int(lines.next().strip())

            # 31
            if (self.scan_mode.upper() == 'REGULAR'):
                block.abscissa_label = lines.next().strip()
                block.abscissa_units = lines.next().strip()
                block.abscissa_start = float(lines.next().strip())
                block.abscissa_increment = float(lines.next().strip())

            # 32
            block.number_of_corresponding_variables = int(lines.next().strip())
            for i in range(block.number_of_corresponding_variables):
                block.corresponding_variable_labels.append(lines.next().strip())
                block.corresponding_variable_units.append(lines.next().strip())

            # 33
            block.signal_mode = lines.next().strip()

            # 34
            block.signal_collection_time = float(lines.next().strip())

            # 35
            block.number_of_scans_to_compile_this_block = int(lines.next().strip())

            # 36
            block.signal_time_correction = float(lines.next().strip())

            # 37
            if ( (block.technique.upper() == 'AES DIFF' or
                  block.technique.upper() == 'AES DIR' or
                  block.technique.upper() == 'EDX' or
                  block.technique.upper() == 'ELS' or
                  block.technique.upper() == 'UPS' or
                  block.technique.upper() == 'XPS' or
                  block.technique.upper() == 'XRF') and
                 (self.experiment_mode.upper() == 'MAPDP' or
                  self.experiment_mode.upper() == 'MAPSVDP' or
                  self.experiment_mode.upper() == 'SDP' or
                  self.experiment_mode.upper() == 'SDPSV') ):
                block.sputtering_source_energy = float(lines.next().strip())
                block.sputtering_source_beam_current = float(lines.next().strip())
                block.sputtering_source_width_x = float(lines.next().strip())
                block.sputtering_source_width_y = float(lines.next().strip())
                block.sputtering_source_polar_angle_of_incidence = float(lines.next().strip())
                block.sputtering_source_azimuth = float(lines.next().strip())
                block.sputtering_mode = lines.next().strip()

            # 38
            block.sample_normal_polar_angle_of_tilt = float(lines.next().strip())
            block.sample_normal_tilt_azimuth = float(lines.next().strip())

            # 39
            block.sample_rotation_angle = float(lines.next().strip())

            # 40
            block.number_of_additional_numerical_parameters = int(lines.next().strip())
            for i in range(block.number_of_additional_numerical_parameters):
                block.additional_numerical_parameter_labels.append(lines.next().strip())
                block.additional_numerical_parameter_units.append(lines.next().strip())
                block.additional_numerical_parameter_values.append(float(lines.next().strip()))

            # end numbered lines in block specification
            for i in range(self.number_of_future_upgrade_block_entries):
                block.future_block_entries.append(lines.next().strip())

            block.number_of_ordinate_values = int(lines.next().strip())

            variables = []
            for i in range(block.number_of_corresponding_variables):
                block.minimum_ordinate_values.append(float(lines.next().strip()))
                block.maximum_ordinate_values.append(float(lines.next().strip()))
                variables.append([])

            for i in range(block.number_of_ordinate_values/block.number_of_corresponding_variables):
                for j in range(block.number_of_corresponding_variables):
                    variables[j].append(float(lines.next().strip()))

            block.ordinate_values = variables

            self.blocks.append(block)
        #end for i in range(self.number_of_blocks):

    #end def read_from_file(self):

    def write_to_file(self):
        self.file.seek(0)
        lines = iter(self.file.readlines())

        self.format_identifier = lines.next().strip()
        self.institution_identifier = lines.next().strip()
        self.instrument_model_identifier = lines.next().strip()
        self.operator_identifier = lines.next().strip()
        self.experiment_identifier = lines.next().strip()
        self.number_of_lines_in_comment = int(lines.next().strip())

        self.comment = ""
        for i in range(self.number_of_lines_in_comment):
            self.comment = self.comment + lines.next()

        self.experiment_mode = lines.next().strip()
        self.scan_mode = lines.next().strip()

        if (self.experiment_mode.upper() == 'MAP' or
            self.experiment_mode.upper() == 'MAPDP' or
            self.experiment_mode.upper() == 'NORM' or
            self.experiment_mode.upper() == 'SDP'):
            self.number_of_spectral_regions = int(lines.next().strip())

        if (self.experiment_mode.upper() == 'MAP' or
            self.experiment_mode.upper() == 'MAPDP'):
            self.number_of_analysis_positions = int(lines.next().strip())
            self.number_of_discrete_x_coordinates_available_in_full_map = int(lines.next().strip())
            self.number_of_discrete_y_coordinates_available_in_full_map = int(lines.next().strip())

        self.number_of_experimental_variables = int(lines.next().strip())

        for i in range(self.number_of_experimental_variables):
            self.experimental_variable_labels.append(lines.next().strip())
            self.experimental_variable_units.append(lines.next().strip())

        self.number_of_entries_in_parameter_inclusion_or_exclusion_list = int(lines.next().strip())

        for i in range(abs(self.number_of_entries_in_parameter_inclusion_or_exclusion_list)):
            self.parameter_inclusion_or_exclusion_prefix_numbers.append(int(lines.next().strip()))

        self.number_of_manually_entered_items_in_block = int(lines.next().strip())

        for i in range(self.number_of_manually_entered_items_in_block):
            self.prefix_numbers_of_manually_entered_items.append(int(lines.next().strip()))

        self.number_of_future_upgrade_experiment_entries = int(lines.next().strip())

        self.number_of_future_upgrade_block_entries = int(lines.next().strip())

        for i in range(self.number_of_future_upgrade_block_entries):
            self.future_upgrade_experiment_entries.append(lines.next().strip())

        self.number_of_blocks = int(lines.next().strip())

        for i in range(self.number_of_blocks):
            block = VAMASBlock()
            block.block_identifier = lines.next().strip()
            block.sample_identifier = lines.next().strip()
            # 1
            year = int(lines.next().strip())

            # 2
            month = int(lines.next().strip())

            # 3
            day = int(lines.next().strip())

            # 4
            hours = int(lines.next().strip())

            # 5
            minutes = int(lines.next().strip())

            # 6
            seconds = int(lines.next().strip())
            block.date = datetime(year,month,day,hours,minutes,seconds)

            # 7
            block.number_of_hours_in_advance_of_greenwich_mean_time = int(lines.next().strip())

            # 8
            block.number_of_lines_in_block_comment = int(lines.next().strip())
            for i in range(block.number_of_lines_in_block_comment):
                block.comment = block.comment + lines.next()

            # 9
            block.technique = lines.next().strip()

            # 10
            if (self.experiment_mode.upper() == 'MAP' or
                self.experiment_mode.upper() == 'MAPDP'):
                block.x_coordinate = int(lines.next().strip())
                block.y_coordinate = int(lines.next().strip())

            # 11
            for i in range(self.number_of_experimental_variables):
                block.values_of_experimental_variables.append(float(lines.next().strip()))

            # 12
            block.analysis_source_label = lines.next().strip()

            # 13
            if (self.experiment_mode.upper() == 'MAPDP' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SDP' or
                self.experiment_mode.upper() == 'SVDP' or
                block.technique.upper() == 'FABMS' or
                block.technique.upper() == 'FABMS energy spec'.upper() or
                block.technique.upper() == 'ISS' or
                block.technique.upper() == 'SIMS' or
                block.technique.upper() == 'SIMS energy spec'.upper() or
                block.technique.upper() == 'SNMS' or
                block.technique.upper() == 'SNMS energy spec'.upper()):
                block.sputtering_ion_or_atomic_number = lines.next().strip()
                block.number_of_atoms_in_sputtering_ion_or_atom_particle = int(lines.next().strip())
                block.sputtering_ion_or_atom_charge_sign_and_number = int(lines.next().strip())

            # 14
            block.analysis_source_characteristic_energy = float(lines.next().strip())

            # 15
            block.analysis_source_strength = float(lines.next().strip())

            # 16
            block.analysis_source_beam_width_x = float(lines.next().strip())
            block.analysis_source_beam_width_y = float(lines.next().strip())

            # 17
            if (self.experiment_mode.upper() == 'MAP' or
                self.experiment_mode.upper() == 'MAPDP' or
                self.experiment_mode.upper() == 'MAPSV' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SEM'):
                block.field_of_view_x = float(lines.next().strip())
                block.field_of_view_y = float(lines.next().strip())

            # 18
            if (self.experiment_mode.upper() == 'MAPSV' or
                self.experiment_mode.upper() == 'MAPSVDP' or
                self.experiment_mode.upper() == 'SEM'):
                block.first_linescan_start_x_coordinate = int(lines.next().strip())
                block.first_linescan_start_y_coordinate = int(lines.next().strip())
                block.first_linescan_finish_x_coordinate = int(lines.next().strip())
                block.first_linescan_finish_y_coordinate = int(lines.next().strip())
                block.last_linescan_finish_x_coordinate = int(lines.next().strip())
                block.last_linescan_finish_y_coordinate = int(lines.next().strip())

            # 19
            block.analysis_source_polar_angle_of_incidence = float(lines.next().strip())

            # 20
            block.analysis_source_azimuth = float(lines.next().strip())

            # 21
            block.analyser_mode = lines.next().strip()

            # 22
            block.analyser_pass_energy_of_retard_ratio_or_mass_resolution = float(lines.next().strip())

            # 23
            if (block.technique.upper() == 'AES DIFF'):
                block.differential_width = float(lines.next().strip())

            # 24
            block.magnification_of_analyser_transfer_lens = float(lines.next().strip())

            # 25
            block.analyser_work_function_or_acceptance_energy_of_atom_or_ion = float(lines.next().strip())

            # 26
            block.target_bias = float(lines.next().strip())

            # 27
            block.analysis_width_x = float(lines.next().strip())
            block.analysis_width_y = float(lines.next().strip())

            # 28
            block.analyser_axis_take_off_polar_angle = float(lines.next().strip())
            block.analyser_axis_take_off_azimuth = float(lines.next().strip())

            # 29
            block.species_label = lines.next().strip()

            # 30
            block.transition_or_charge_state_label = lines.next().strip()
            block.charge_of_detected_particle = int(lines.next().strip())

            # 31
            if (self.scan_mode.upper() == 'REGULAR'):
                block.abscissa_label = lines.next().strip()
                block.abscissa_units = lines.next().strip()
                block.abscissa_start = float(lines.next().strip())
                block.abscissa_increment = float(lines.next().strip())

            # 32
            block.number_of_corresponding_variables = int(lines.next().strip())
            for i in range(block.number_of_corresponding_variables):
                block.corresponding_variable_labels.append(lines.next().strip())
                block.corresponding_variable_units.append(lines.next().strip())

            # 33
            block.signal_mode = lines.next().strip()

            # 34
            block.signal_collection_time = float(lines.next().strip())

            # 35
            block.number_of_scans_to_compile_this_block = int(lines.next().strip())

            # 36
            block.signal_time_correction = float(lines.next().strip())

            # 37
            if ( (block.technique.upper() == 'AES DIFF' or
                  block.technique.upper() == 'AES DIR' or
                  block.technique.upper() == 'EDX' or
                  block.technique.upper() == 'ELS' or
                  block.technique.upper() == 'UPS' or
                  block.technique.upper() == 'XPS' or
                  block.technique.upper() == 'XRF') and
                 (self.experiment_mode.upper() == 'MAPDP' or
                  self.experiment_mode.upper() == 'MAPSVDP' or
                  self.experiment_mode.upper() == 'SDP' or
                  self.experiment_mode.upper() == 'SDPSV') ):
                block.sputtering_source_energy = float(lines.next().strip())
                block.sputtering_source_beam_current = float(lines.next().strip())
                block.sputtering_source_width_x = float(lines.next().strip())
                block.sputtering_source_width_y = float(lines.next().strip())
                block.sputtering_source_polar_angle_of_incidence = float(lines.next().strip())
                block.sputtering_source_azimuth = float(lines.next().strip())
                block.sputtering_mode = lines.next().strip()

            # 38
            block.sample_normal_polar_angle_of_tilt = float(lines.next().strip())
            block.sample_normal_tilt_azimuth = float(lines.next().strip())

            # 39
            block.sample_rotation_angle = float(lines.next().strip())

            # 40
            block.number_of_additional_numerical_parameters = int(lines.next().strip())
            for i in range(block.number_of_additional_numerical_parameters):
                block.additional_numerical_parameter_labels.append(lines.next().strip())
                block.additional_numerical_parameter_units.append(lines.next().strip())
                block.additional_numerical_parameter_values.append(float(lines.next().strip()))

            # end numbered lines in block specification
            for i in range(self.number_of_future_upgrade_block_entries):
                block.future_block_entries.append(lines.next().strip())

            block.number_of_ordinate_values = int(lines.next().strip())

            variables = []
            for i in range(block.number_of_corresponding_variables):
                block.minimum_ordinate_values.append(float(lines.next().strip()))
                block.maximum_ordinate_values.append(float(lines.next().strip()))
                variables.append([])

            for i in range(block.number_of_ordinate_values/block.number_of_corresponding_variables):
                for j in range(block.number_of_corresponding_variables):
                    variables[j].append(float(lines.next().strip()))

            block.ordinate_values = variables

            self.blocks.append(block)

        #end for i in range(self.number_of_blocks):
    #end write_to_file
#end Class VAMASExperiment

#Class to store the data contained within one block of an experiment
class VAMASBlock:
    def __init__(self):
        self.set_defaults()

    def abscissa(self):
        num_values = self.number_of_ordinate_values/self.number_of_corresponding_variables
        if (self.technique.upper() != 'AES DIR'):
            if (self.abscissa_label.upper() == 'BINDING ENERGY'):
                return -(arange(num_values) * self.abscissa_increment + self.abscissa_start)
            else:
                return -(self.analysis_source_characteristic_energy - (arange(num_values) * self.abscissa_increment + self.abscissa_start))
        else:
            return arange(num_values) * self.abscissa_increment + self.abscissa_start

    def ordinate(self,variable_number):
        return array(self.ordinate_values[variable_number])

    def set_defaults(self):
        self.block_identifier = "Not Specified"
        self.sample_identifier = "Not Specified"
        self.date = datetime.now()
        self.number_of_hours_in_advance_of_greenwich_mean_time = 0
        self.number_of_lines_in_block_comment = 0
        self.comment = "Not Specified"
        self.technique = "XPS"
        # ( 'AES diff | 'AES dir' | 'EDX' | 'ELS' | 'FABMS' | 'FABMSenergy spec' | 'ISS' | 'SIMS' |
        #   'SIMS energy spec' | 'SNMS' | 'SNMS energy spec' | 'UPS' | 'XPS' | 'XRF' )

        # The ordinal numbers, starting with unity, of the point in the array along the analysis source deflection system x,y-axes.
        # The following two entries are inserted if and only if the value of experiment mode is either 'MAP' or 'MAPDI''.
        self.x_coordinate = 0
        self.y_coordinate = 0

        # value of experimental variable may be, for example, total time in seconds,
        # total etch time in seconds, temperature in Kelvin, energy in electron volts
        # or mass in unified atomic mass units.
        # The number of occurrences of value of experimental variable is specified by
        # the value of number of experimental variables above, and the order in which
        # the values are given is the same as the order in which experimental variable
        # label and experimental variable units are declared above.
        self.values_of_experimental_variables = []

        self.analysis_source_label = "Not Specified"

        # The following three entries are inserted if and only if either
        # (1) the value of experiment mode is 'MAPDP, 'MAPSVDP, 'SDP or 'SDPSV, or
        # (2) the value of technique is 'FABMS', 'FABMS energy spec', 'ISS', 'SIMS',
        #     'SIMSenergy spec', 'SNMS' or 'SNMS energy spec'.
        self.sputtering_ion_or_atomic_number = 0
        self.number_of_atoms_in_sputtering_ion_or_atom_particle = 0
        self.sputtering_ion_of_atom_charge_sign_and_number = -1

        # energy in electron volts
        self.analysis_source_characteristic_energy = 0

        # power in watts for XPS and XRF; beam current in nanoamps for AES, EDX, ISS, SIMS, and SNMS; beam equivalent for FABMS
        self.analysis_source_strength = 0

        # width in micrometres at the sample in the plane perpendicular to the source beam (x,y)
        self.analysis_source_beam_width_x = 0
        self.analysis_source_beam_width_y = 0

        # FOV in micrometers
        # These two entries are inserted if and only if the value of experiment mode is 'MAP', 'MAPDP', 'MAPSV', 'MAPSVDP' or 'SEM'.
        self.field_of_view_x = 0
        self.field_of_view_y = 0

        # The followwing six entries are inserted if and only if the value of experiment mode
        # is 'MAPSV', 'MAPSVDP' or 'SEM'.
        # They are required for specifymg the size and shape of the map and for relating the
        # order in the scan sequence to the position on the sample.
        # In the coordinate system to be used, x-values start at unity at the left-hand side
        # of the frame and increase towards the right-hand side, and y-values start at unity
        # at the top of the frame and increase towards the bottom of the frame, as shown below.
        #
        #  +-------------------+
        #  |1,1             N,1|
        #  |                   |
        #  |                   |
        #  |                   |
        #  |1,M             N,M|
        #  +-------------------+
        #
        self.first_linescan_start_x_coordinate = 0
        self.first_linescan_start_y_coordinate = 0
        self.first_linescan_finish_x_coordinate = 0
        self.first_linescan_finish_y_coordinate = 0
        self.last_linescan_finish_x_coordinate = 0
        self.last_linescan_finish_y_coordinate = 0

        # degrees from upward zdirection, defined by the sample stage
        self.analysis_source_polar_angle_of_incidence = 0

        # degrees clockwise from the y-direction towards the operator, defined by the sample stage
        self.analysis_source_azimuth = 0

        self.analyser_mode = "FAT"
        # ( 'FAT | 'FRR' | 'constant delta m' | 'constant m/delta m' )

        # energy in electron volts, mass in amu
        self.analyser_pass_energy_or_retard_ratio_or_mass_resolution = 0

        # electron volts peak-to-peak for sinusoidal modulation or computer differentiation.
        # differential width is inserted if and only if the value of technique is 'AES diff'.
        self.differential_width = 0

        self.magnification_of_analyser_transfer_lens = 1

        # positive value for work function in electron volts for AES,ELS, ISS, UPS and XPS.
        # The acceptance energy of an ion is the energy filter pass energy of the mass spectrometer for FABMS, SIMS, and SNMS.
        self.analyser_work_function_or_acceptance_energy_of_atom_or_ion = 0

        # target bias is in volts, including the sign
        self.target_bias = 0

        # The analysis width x is the gated signal width of the source in the x-direction in the
        # plane perpendicular to the beam for FABMS, FABMS energy spec, ISS, SIMS, SIMS energy spec,
        # SNMS and SNMS energy spec, the analyser slit length divided by the magnification of the
        # analyser transfer lens to that slit for AES diff, AES dir, ELS, UPS and XPS, and is the
        # source width in the x-direction for both EDX and XRF.
        #  analysis width is in micrometres.
        self.analysis_width_x = 0
        self.analysis_width_y = 0

        # degrees from upward z-direction, defined by the sample stage
        self.analysis_axis_take_off_polar_angle = 0

        # degrees clockwise from the y-direction towards the operator, defined by the sample stage
        self.analyser_axis_take_off_azimuth = 0

        # elemental symbol or molecular formula
        self.species_label = "Not Specified"

        # example: 'KLL' for AES, '1s' for XPS, '-1' for SIMS
        self.transition_or_charge_state_label = "Not Specified"

        # example: -1 for AES and XPS, +1 for positive SIMS
        self.charge_of_detected_particle = -1

        # the following four entries are inserted if and only if the value of scan mode is 'REGULAR'.
        self.abscissa_label = "Not Specified"  # text line
        self.abscissa_units = "Not Specified"  # units
        self.abscissa_start = 0      # real number
        self.abscissa_increment = 0  # real number

        # If the data is in the form of sets of corresponding values of two or more variables then
        # number_of_corresponding_variables is equal to the number of variables, otherwise it is equal to unity.
        self.number_of_corresponding_variables = 1

        # The number of occurrences of the below pair of items is specified by the value of
        # number_of_corresponding_variables above.
        self.corresponding_variable_labels = []
        self.corresponding_variable_units = []

        self.signal_mode = "Not Specified"
        # ( 'analogue' | 'pulse counting' ),
        # Analogue signals, while recorded digitally, may be of either sign and have a gain
        # which may be noted in the block comment. Pulse counting signals are integers with
        # values equal to or greater than zero.


        # time in seconds per scan for each channel or array-point, except for both EDX and XRF
        # where it is the total spectrum collection time
        self.signal_collection_time = 0

        self.number_of_scans_to_compile_this_block = 1

        # This is the system dead time, except for EDX and XRF where it is the livetime-corrected
        # acquisition time. In the case of a dead time, a positive value indicates that the count
        # rate should be corrected by dividing by (1 - measured rate x dead time) whereas a negative
        # value indicates a correction by multiplying by (exp(1rue count rate x dead time)).
        #    signal_time_correction is in seconds
        self.signal_time_correction = 0

        # The following seven entries are for a sputtering source used in addition to the
        # analysis source, as in depth profiling, in AES diff, AES dir, EDX, ELS, UPS, XPS or XRF.
        # The following seven entries are inserted if and only if both
        #  (1) the value of technique is 'AES diff, 'AES dir', 'EDX, 'ELS', 'UPS', 'XPS' or 'XRF, and
        #  (2) the value of experiment mode is 'MAPDP, 'MAPSVDP, 'SDP or 'SDPSV'

        # energy in electron volts
        self.sputtering_source_energy = 0

        # current in nanoamps or equivalent for neutrals
        self.sputtering_source_beam_crrent = 0

        # width in micrometres at the sample in the plane perpendicular to the sputtering source beam x,y
        self.sputtering_source_width_x = 0
        self.sputtering_source_width_y = 0

        #degrees from upward z-direction, defined by the sample stage
        self.sputtering_source_polar_angle_of_incidence = 0

        # degrees clockwise from the y-direction towards the operator, defined by the sample stage
        self.sputtering_source_polar_angle_of_incidence = 0

        # The value of sputtering mode is either 'continuous', when sputtering continues while
        # spectral data is being recorded, or 'cyclic', when sputtering is suspended while spectral data is being recorded.
        self.sputtering_mode = "Not Specified"
        # ( 'continuous' | 'cyclic' )

        # degrees from upward z-direction, defined by the sample stage
        self.sample_normal_polar_angle_of_tilt = 0

        # degrees clockwise from the y-direction towards the operator, defined by the sample stage
        self.sample_normal_tilt_azimuth = 0

        # degrees clockwise rotation about the sample normal. If this is referenced to a
        # particular direction on the sample this direction would be specified in a comment line at item number 8.
        self.sample_rotation_angle = 0

        self.number_of_additional_numerical_parameters = 0

        # The number of occurrences of the following group of three entries is specified by the value of
        # number_of_additional_numerical_parameters.
        self.additional_numerical_parameter_labels = []   # text line
        self.additional_numerical_parameter_units = []    # units
        self.additional_numerical_parameter_values = []   # real number

        # The number of occurrences of future upgrade block entry is given by the value of number of
        # future upgrade block entries above. It is defined as a text line so that any integer, real
        # number or text line inserted at this point by a future upgrade of the Format can be read
        # as a text line then discarded.
        self.future_upgrade_block_entries = []

        # The value of number of ordinate values is equal to product of the value of number of
        # corresponding variables and the number of sets of corresponding variables to be transferred.
        self.number_of_ordinate_values = 1

        # The number of occurrences of the following pair of entries is specified by the value of
        # number_of_corresponding_variables above. The order in which the pairs of entries
        # appear is the same as the order in which the corresponding values of
        # corresponding_variable_label are given above.
        self.minimum_ordinate_values = []
        self.maximum_ordinate_values = []

        # The number of occurrences of ordinate value is specified by the value of number of
        # ordinate values above. If the value of number of corresponding variables is greater
        # than unity then the data is sent in the form of successive complete sets, each set
        # consisting of an ordinate value for each of the corresponding variables arranged in
        # the same order as that in which each value of corresponding variable label is given above.
        #   There must be at least one ordinate value.
        self.ordinate_values = []
