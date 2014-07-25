from pylab import *

def read_bl7011_xas(filename):

    experiment = np.genfromtxt(filename,
        dtype='f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8',
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

    return experiment

def read_augerscan(filename):

    experiment = np.genfromtxt(filename,
        dtype='f8,f8',
        names=['Energy', 'Counts'],
        skip_header=1)

    return experiment
