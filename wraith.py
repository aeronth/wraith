#!/usr/bin/env ipython

import sys

from PySide.QtCore import QCoreApplication
from PySide.QtGui import QApplication, QStyleFactory

import matplotlib
matplotlib.rcParams['backend.qt4'] = 'PySide'

from wraith_main_window import Form

#main function to start up program
def main():
    matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
    app = QApplication(sys.argv)
    form = Form()
    #style choices common to both modes        
    QApplication.setStyle(QStyleFactory.create('Plastique'))
    QApplication.setPalette(QApplication.style().standardPalette())
    #do it
    form.show()
    app.exec_()

#wraith function to start up program from interactive terminal
def wraith():
    matplotlib.rcParams['mathtext.fontset'] = 'stixsans'
    app = QCoreApplication.instance()
    app.form = Form()
    #style choices common to both modes    
    QApplication.setStyle(QStyleFactory.create('Plastique'))
    QApplication.setPalette(QApplication.style().standardPalette())
    #do it
    app.form.show()

#if run from commandline then start up by calling main()
if __name__ == "__main__":
    main()
else:
    app = QCoreApplication.instance()
