'''
Created on 18-03-2012

@author: Pawel Jarosz
'''
import os, sys
import atexit

from PySide import QtCore, QtGui

from IPython.zmq.ipkernel import IPKernelApp
from IPython.lib.kernel import find_connection_file, connect_qtconsole
from IPython.frontend.qt.kernelmanager import QtKernelManager
from IPython.frontend.qt.console.rich_ipython_widget import RichIPythonWidget
from IPython.config.application import catch_config_error

class IPythonLocalKernelApp(IPKernelApp):
    """IPython kernel application with nonblocking loop, running in dedicated thread.
    example:
        app = QtGui.QApplication([])
        kernelapp = IPythonLocalKernelApp.instance()
        kernelapp.start()
        namespace = kernelapp.get_user_namespace()
        namespace["QtGui"]=QtGui
        namespace["QtCore"]=QtCore
        app.exec_()"""
    #DEFAULT_INSTANCE_ARGS starting commandline
    DEFAULT_INSTANCE_ARGS = ['qtconsole','--pylab=inline', '--colors=linux']

    @catch_config_error
    def initialize(self, argv=None):
        super(IPythonLocalKernelApp, self).initialize(argv)
        self.kernel.eventloop = self.loop_qt4_nonblocking
