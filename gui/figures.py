from PyQt4 import QtGui, QtCore, uic
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt


def create_histogram():
	canvas = HistoCanvas()
	return canvas


class HistoCanvas(FigureCanvas):

    def __init__(self, dpi=100):

        self.fig = Figure(dpi=dpi)                
        self.axes = self.fig.add_subplot(111)        
        FigureCanvas.__init__(self, self.fig)
        FigureCanvas.setSizePolicy(self,
                                   QtGui.QSizePolicy.Expanding,
                                   QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

	def drawGraph(self, data):
		self.axes.cla()
		self.axes.hist(data)
		self.draw()