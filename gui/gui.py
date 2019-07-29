import pickle, os
from PyQt4 import QtGui, QtCore, uic
import sys
sys.path.insert(0, "..")
from procell import Simulator
from numpy import loadtxt
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from numpy.random import uniform
from pylab import *
import pandas as pd
from PyQt4.QtCore import QThread, pyqtSignal
import resources
from estimator import hellinger1
from project import Project
from PyQt4.QtGui import QApplication
try:
	import seaborn as sns  
except:
	print "WARNING: Seaborn not installed"
from estimator import Calibrator

def dummy(): return 0

matplotlib.rcParams.update({'font.size': 8})

def spawn():
	from subprocess import check_output
	new_exec = [sys.executable, os.path.realpath(__file__), "empty"]
	ret = check_output(new_exec)

def rebin(series, lower, upper, N=1000):
	bins = logspace(ceil(log10(lower)), ceil(log10(upper)), N)
	res  = zeros(N, dtype=int)
	pos  = 0
	for k,v in zip(series.T[0], series.T[1]):
		while(k>bins[pos]):
			pos+=1
		res[pos] += v
	return res, bins

def verticalResizeTableViewToContents(tableView):
    count=tableView.verticalHeader().count()
    scrollBarHeight=tableView.horizontalScrollBar().height()
    horizontalHeaderHeight=tableView.horizontalHeader().height()
    rowTotalHeight = 0
    for i in xrange(count):
        rowTotalHeight+=tableView.verticalHeader().sectionSize(i)
    tableView.setMinimumHeight(horizontalHeaderHeight+rowTotalHeight+scrollBarHeight)


class PandasModel(QtCore.QAbstractTableModel):
    """
    Class to populate a table view with a pandas dataframe
    """
    def __init__(self, data, parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent)
        self._data = data

    def rowCount(self, parent=None):
        return len(self._data.values)

    def columnCount(self, parent=None):
        return self._data.columns.size

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if index.isValid():
            if role == QtCore.Qt.DisplayRole or role == QtCore.Qt.EditRole:
                return str(self._data.values[index.row()][index.column()])
        return None

    def headerData(self, col, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self._data.columns[col]
        return None

 	def sort(self, column, order):
		colname = self._df.columns.tolist()[column]
		self.layoutAboutToBeChanged.emit()
		self._df.sort_values(colname, ascending= order == QtCore.Qt.AscendingOrder, inplace=True)
		self._df.reset_index(inplace=True, drop=True)
		self.layoutChanged.emit()

	def setData(self, index, value, role=QtCore.Qt.EditRole):
		if index.isValid():
			self._data.iloc[index.row(), index.column()] = value
			if self.data(index, QtCore.Qt.DisplayRole) == value:
				self.dataChanged.emit(index, index)
				return True
		return False

	def flags(self, index):
		return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable | QtCore.Qt.ItemIsEditable

class AboutWindow(QtGui.QDialog):

	def __init__(self):
		super(AboutWindow, self).__init__()
		uic.loadUi('about.ui', self)
		


class MainWindow(QtGui.QMainWindow):

	def __init__(self):
		super(MainWindow, self).__init__()
		uic.loadUi('mainwindow.ui', self)

		self._about_window = AboutWindow()

		self._initial_histo    = None
		self._target_histo     = None
		self._validation_histo = None
		self._simulated_histo  = None
		self._simulated_validation_histo = None

		self._initial_histo_path = None
		self._target_histo_path = None
		self._validation_histo_path = None

		self._initial_histo_figure    = plt.figure(figsize=(10,7), tight_layout=True)
		self._target_histo_figure     = plt.figure(figsize=(10,7), tight_layout=True)
		self._validation_histo_figure = plt.figure(figsize=(10,7), tight_layout=True)

		self._initial_histo_ax    = self._initial_histo_figure.add_subplot(111)
		self._target_histo_ax     = self._target_histo_figure.add_subplot(111)
		self._validation_histo_ax = self._validation_histo_figure.add_subplot(111)

		self._initial_histo_canvas    = FigureCanvas(self._initial_histo_figure)
		self._target_histo_canvas     = FigureCanvas(self._target_histo_figure)
		self._validation_histo_canvas = FigureCanvas(self._validation_histo_figure)
	
		self._initial_histo_canvas.setStyleSheet("background:transparent;"); self._initial_histo_figure.set_facecolor("#ffff0000")
		self._target_histo_canvas.setStyleSheet("background:transparent;"); self._target_histo_figure.set_facecolor("#ffff0000")
		self._validation_histo_canvas.setStyleSheet("background:transparent;"); self._validation_histo_figure.set_facecolor("#ffff0000")

		self._target_histo_ax.set_xlabel("Fluorescence")
		self._target_histo_ax.set_ylabel("Cells frequency")

		self._initial_histo_ax.set_xlabel("Fluorescence")
		self._initial_histo_ax.set_ylabel("Cells frequency")

		self.initial_layout.addWidget (self._initial_histo_canvas)
		self.target_layout.addWidget (self._target_histo_canvas)
		self.validation_layout.addWidget (self._validation_histo_canvas)

		self._population_names 			= []
		self._population_proportions 	= []
		self._population_means 			= []
		self._population_std   			= []
		self._population_info  			= []
		self._population_minmean 		= []
		self._population_maxmean 		= []
		self._population_minsd 			= []
		self._population_maxsd 			= []

		# debug
		"""
		if len(sys.argv)<2:
			self._add_population("quiescent", proportion=0.06, mean="-", st="-")
			self._add_population("slowly proliferating", 0.47, mean=63.7, st=28.9)
			self._add_population("fast proliferating", 0.47, mean=41.1, st=14.3)
			self._import_initial_histo("C:\\AML9_Exp3_t10d_UT1_BM_GFPpos.txt")
			self._import_target_histo("C:\\AML9_Exp3_t10d_DOX3_BM_GFPpos.txt")
			self._import_validation_histo("C:\\AML9_Exp3_t21d_DOX_BM_GFPpos.txt")
		"""


		self._update_populations()
		self.populations_table.setColumnWidth(0, 160) # name
		self.populations_table.setColumnWidth(1, 120) # prop
		self.populations_table.setColumnWidth(2, 130) # mu
		self.populations_table.setColumnWidth(3, 130) # sigma
		self.populations_table.setColumnWidth(4, 200) # minmean
		self.populations_table.setColumnWidth(5, 200) # maxmean
		self.populations_table.setColumnWidth(6, 270) # minsd
		self.populations_table.setColumnWidth(7, 270) # maxsd
		self.populations_table.setColumnWidth(8, 160) # info
		verticalResizeTableViewToContents(self.populations_table)

		self.Simulator = Simulator()

		# status bar
		self.statusBar = QtGui.QStatusBar()
		self.setStatusBar(self.statusBar)
		self.statusBar.showMessage("ProCell ready.")

		# progress bar
		self.progress = QtGui.QProgressBar()
		self.progress.setRange(0,100)
		self.statusBar.addPermanentWidget(self.progress)

		"""
		# project stuff
		#self._load_project_from_file("./prova.prc")		
		"""
		self._project_filename = None
		self._version = " 1.2.0"

		self._update_window_title()

		self._calibrator = Calibrator()
		self._calibrator.distribution = "gauss"  # default semantics

		self._load_project_from_file("../models/model3.prc")

		# implement in preferences (TODO)
		self._path_to_GPU_procell = "D:\\GITHUBrepos\\cuda-pro-cell\\build\\bin\\Release\\procell.exe"

		self.show()

	def _new_project(self):
		spawn()


	def _show_about(self):
		self._about_window.show()

	def _update_populations(self):
		data = self._get_population_data_frame()
		self.model = PandasModel(data=data)	
		self.populations_table.setModel(self.model)
		

	def _add_population(self, name, proportion, mean, st, minimum_mean=0, maximum_mean=0, minimum_sd=0, maximum_sd=0, info=None):
		self._population_names.append(name)
		self._population_proportions.append(proportion)
		self._population_means.append(mean)
		self._population_std.append(st)
		self._population_info.append(info)
		self._population_minmean.append(minimum_mean)
		self._population_maxmean.append(maximum_mean)
		self._population_minsd.append(minimum_sd)
		self._population_maxsd.append(maximum_sd)

		verticalResizeTableViewToContents(self.populations_table)

	def _get_population_data_frame(self):
		return  pd.DataFrame({	'1. Population name': self._population_names, 
								'2. Proportion': self._population_proportions, 
								'3. Mean division time': self._population_means, 
								'4. Standard deviation': self._population_std, 
								'5. Mean division time (minimum)': self._population_minmean,
								'6. Mean division time (maximum)': self._population_maxmean,
								'7. Standard deviation division time (minimum)': self._population_minsd,
								'8. Standard deviation division time (maximum)': self._population_maxsd,
								'9. Info': self._population_info}
						)

	def drop_histogram(self):
		print "NOT IMPLEMENTED YET"

	def _message_error(self, text, additional=None):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Critical)
		msg.setText(text)
		msg.setWindowTitle("Import error")
		if additional!=None:
			msg.setInformativeText("This is additional information")
			msg.setDetailedText(additional)
		msg.setStandardButtons(QtGui.QMessageBox.Ok)
		msg.exec_()

	def _import_initial_histo(self, path):
		path = str(path)
		try:
			A = loadtxt(path)
			#print "Valid histogram found. Data:"; 		print A
			self._initial_histo = A[:]
			self._initial_histo_path = path
			self._update_initial_plot()		
			return True

		except ValueError:
			self._message_error("ERROR: invalid histogram specified"); 
			self._initial_histo = None
			return False

	def _import_target_histo(self, path):
		path = str(path)
		try:
			A = loadtxt(path)
			#print "Valid histogram found. Data:"; 			print A
			self._target_histo = A[:]
			self._target_histo_path = path
			self._update_target_plot()
			self.normtotarget.setEnabled(True)
			return True

		except ValueError:
			self._message_error("ERROR: invalid histogram specified")
			self._target_histo = None
			self.normtotarget.setEnabled(False)
			return False

	def _import_validation_histo(self, path):
		path = str(path)
		try:
			A = loadtxt(path)
			#print "Valid histogram found. Data:"; 			print A
			self._validation_histo = A[:]
			self._validation_histo_path = path
			self._update_validation_plot()		
			return True

		except ValueError:
			self._message_error("ERROR: invalid histogram specified")
			self._validation_histo = None
			return False

	def import_target_histo(self):
		# open dialog
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Open histogram file', 'c:\\',"Histogram files (*.fcs *.txt)")

		# import file
		res = self._import_target_histo(fname)

	def import_validation_histo(self):
		# open dialog
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Open histogram file', 'c:\\',"Histogram files (*.fcs *.txt)")

		# import file
		res = self._import_validation_histo(fname)

	def _update_all_plots(self):
		if self._initial_histo 		is not None:	self._update_initial_plot()
		if self._target_histo  		is not None: self._update_target_plot()
		if self._validation_histo 	is not None: self._update_validation_plot()

	def _update_target_plot(self):
		self._target_histo_ax.cla()
		#self._target_histo_ax = self._target_histo_figure.add_subplot(111)

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())

		if self.usecolors.isChecked():
			selected_color = "cyan"
		else:
			selected_color = "gray"

		if self._target_histo is not None:
			res, bins = rebin(self._target_histo, lower, higher, N=calcbins)
			self._target_histo_ax.bar(bins[1:-1], res[1:-1], width=diff(bins[1:]),  color=selected_color, label="Target", alpha=0.5, ec="black", linewidth=0.4, align="edge")
			#print " * Integral of target:", sum(res)

		if self.usecolors.isChecked():
			selected_color = "green"
		else:
			selected_color = "black"

		if self._simulated_histo is not None:
			res2, bins2 = rebin(self._simulated_histo, lower, higher, N=calcbins)
			if self.normtotarget.isChecked():
				ratio = 1.*sum(res2)/sum(res)
			else:
				ratio = 1.0

			self.hellingertarget.setText( "%.3f"  % hellinger1(res2/ratio, res) )
			self._target_histo_ax.bar(bins2[1:-1], res2[1:-1]/ratio, width=diff(bins[1:]),  color=selected_color, label="Simulation", alpha=0.5, ec="black", linewidth=0.4, align="edge")

		self._target_histo_ax.set_xscale("symlog")
		self._target_histo_ax.set_xlim(bins[0],bins[-1])
		self._target_histo_figure.legend(framealpha=1.)
		self._target_histo_ax.set_xlabel("Fluorescence")
		self._target_histo_ax.set_ylabel("Cells frequency")
		self._target_histo_canvas.draw()


	def import_initial_histo(self):
		# open dialog
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Open histogram file', 'c:\\',"Histogram files (*.fcs *.txt)")

		# import file
		if fname!="":
			self._import_initial_histo(fname)
			
	def _update_initial_plot(self):
		self._initial_histo_ax.clear()

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())
		
		res, bins = rebin(self._initial_histo, lower, higher, N=calcbins)

		self._initial_histo_ax.set_xlabel("Fluorescence")
		self._initial_histo_ax.set_ylabel("Cells frequency")
		if self.usecolors.isChecked():
			selected_color = "red"
		else:
			selected_color = "gray"
		self._initial_histo_ax.bar(bins[1:-1], res[1:-1], width=diff(bins[1:]),  
				color=selected_color, label="Initial histogram", alpha=0.5, linewidth=0.4, ec="black", align="edge")
		
		self._initial_histo_ax.set_xscale("symlog")
		self._initial_histo_ax.set_xlim(bins[0],bins[-1])
		self._initial_histo_figure.legend(framealpha=1.)
		self._initial_histo_canvas.draw()

	def _update_validation_plot(self):
		self._validation_histo_ax.clear()

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())
		
		res, bins = rebin(self._validation_histo, lower, higher, N=calcbins)

		self._validation_histo_ax.set_xlabel("Fluorescence")
		self._validation_histo_ax.set_ylabel("Cells frequency")

		if self.usecolors.isChecked():
			selected_color = "purple"
		else:
			selected_color = "gray"
				
		self._validation_histo_ax.bar(bins[1:-1], res[1:-1], width=diff(bins[1:]),  color=selected_color, label="Validation histogram", alpha=0.5, linewidth=0.4, ec="black", align="edge")

		if self._simulated_validation_histo is not None:
			res2, bins2 = rebin(self._simulated_validation_histo, lower, higher, N=calcbins)
			if self.normtotarget.isChecked():
				ratio = 1.*sum(res2)/sum(res)
			else:
				ratio = 1.0

			self.hellingervalidation.setText( "%.3f"  % hellinger1(res2/ratio, res) )
			self._validation_histo_ax.bar(bins2[1:-1], res2[1:-1]/ratio, width=diff(bins[1:]),  color="green", label="Simulation", alpha=0.5, ec="black", linewidth=0.4, align="edge")

		self._validation_histo_ax.set_xscale("symlog")
		self._validation_histo_ax.set_xlim(bins[0],bins[-1])
		self._validation_histo_figure.legend(framealpha=1.)
		self._validation_histo_canvas.draw()

	def new_population(self):
		self._add_population("unnamed population", proportion=0, mean=0,  st=0, minimum_mean=0, maximum_mean=0, minimum_sd=0, maximum_sd=0, info="")
		self._update_populations()	

	def remove_population(self):
		index = self.populations_table.selectionModel().selectedRows()[0]
		row = index.row()
		print " * Removing row %d..." % row
		del self._population_names[row]
		del self._population_proportions[row]
		del self._population_means[row]
		del self._population_std[row]
		del self._population_minmean[row]
		del self._population_maxmean[row]
		del self._population_minsd[row]
		del self._population_maxsd[row]
		del self._population_info[row]
		self._update_populations()
		verticalResizeTableViewToContents(self.populations_table)


	def _create_new_population_interface(self):
		return QtGui.QWidget()

	def save_figures(self):
		file = str(QtGui.QFileDialog.getExistingDirectory(self, "Select Directory"))
		if file is not None:
			self._initial_histo_figure.savefig(file+os.sep+"initial_histogram.pdf", figsize=(10,7), dpi=300)
			self._target_histo_figure.savefig(file+os.sep+"target_histogram.pdf", figsize=(10,7), dpi=300)
			self._validation_histo_figure.savefig(file+os.sep+"validation_histogram.pdf", figsize=(10,7), dpi=300)

	def _ready_to_simulate(self):
		if self._initial_histo is not None:
			if len(self._population_names)>0:
				return True
			else:
				print "WARNING: no populations specified, aborting."
		else:
			print "WARNING: initial histogram not loaded, aborting."
		return False

	def _ready_to_optimize(self):
		if self._target_histo is not None:
			if len(self._population_names)>0:
				return True
			else:
				print "WARNING: no populations specified, aborting."
		else:
			print "WARNING: target histogram not loaded, aborting."
		return False

	def _check_proportions(self):
		return  (sum(self._population_proportions)-1.0)<1e-3


	def run_simulation(self):
		if self._ready_to_simulate():
			if self._check_proportions():

				self.YTN = SimulationThread(self)
				self.YTN._what="target"
				
				self.connect(self.YTN, QtCore.SIGNAL("finished()"), self._done_simulation)

				self.launch_simulation.setEnabled(False)
				self.runvalidation.setEnabled(False)
				self.abort.setEnabled(True)
				self.YTN.start()
			else:
				self._error_proportions()
		else:
			self._query_for_initial()

	def _error_proportions(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Critical)
		msg.setText("The sum of proportions is not 1. Please check the proportion of cells.")
		msg.setWindowTitle("Unable to run simulation")
		#msg.setStandardButtons(QtGui.QMessageBox.OK)
		ret = msg.exec_()

	def _query_for_initial(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Information)
		msg.setText("Cannot simulate without an initial fluorescence histogram.\nDo you want to import an histogram now?")
		msg.setWindowTitle("Unable to run simulation")
		msg.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.Cancel)
		ret = msg.exec_()
		if ret == QtGui.QMessageBox.Yes:
			self.import_initial_histo()
			return True
		return False

	def _query_for_target(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Information)
		msg.setText("Cannot fit parameters without a target fluorescence histogram.\nDo you want to import an histogram now?")
		msg.setWindowTitle("Unable to run optimization")
		msg.setStandardButtons(QtGui.QMessageBox.Yes | QtGui.QMessageBox.Cancel)
		ret = msg.exec_()
		if ret == QtGui.QMessageBox.Yes:
			self.import_target_histo()
			return True
		return False
			
	def run_validation(self):
		if self._ready_to_simulate():
			self.YTN = SimulationThread(self)
			self.YTN._what="validation"
			
			self.connect(self.YTN, QtCore.SIGNAL("finished()"), self._done_simulation)

			self.launch_simulation.setEnabled(False)
			self.runvalidation.setEnabled(False)
			self.abort.setEnabled(True)
			self.YTN.start()

	def _done_simulation(self):
		result_simulation = self.YTN.result_simulation

		if not self.Simulator._abort_variable:

			sorted_res = array(sorted([ [a,b] for (a,b) in zip(result_simulation.keys(), result_simulation.values())]))

			if self.YTN._what == "target":
				self._simulated_histo = sorted_res
				self.statusBar.showMessage("Simulation completed.")
				self._update_target_plot()
			else:
				self._simulated_validation_histo = sorted_res
				self.statusBar.showMessage("Validation completed.")
				self._update_validation_plot()
		else:
			self.statusBar.showMessage("Simulation aborted.")

		self.launch_simulation.setEnabled(True)
		self.runvalidation.setEnabled(True)
		self.abort.setEnabled(False)
		self.Simulator._abort_variable = False
		self.YTN.timer.stop()
		self.progress.reset()

		#del self.Simulator

	
	def _cancel_simulation(self):
		self.YTN.stop()

	def edit_cell(self, index):
		column_name = self.populations_table.model().headerData( index.column(), QtCore.Qt.Horizontal, QtCore.Qt.DisplayRole  )
		current_value = self.model.data(index)

		# population name
		if index.column()==0: 
			item, ok  = QtGui.QInputDialog.getText(self, "Change population name", "Enter new name:", QtGui.QLineEdit.Normal, current_value )
			if ok and item:
				self._population_names[index.row()] = item
				self._update_populations()

		# proportion
		elif index.column()==1:
			item, ok  = QtGui.QInputDialog.getDouble(self, "Change proportion of population", "Enter new proportion:", float(current_value), 0, 1,3)
			if ok and item:
				self._population_proportions[index.row()] = item
				self._update_populations()

		# mean division time
		elif index.column()==2:
			# item, ok  = QtGui.QInputDialog.getDouble(self, "Change division time", "Enter new mean division time:", float(current_value), 0)
			item, ok  = QtGui.QInputDialog.getText(self, "Change division time", "Enter new mean division time:", QtGui.QLineEdit.Normal, current_value)
			if ok and item:
				try:
					value = float(item)
					self._population_means[index.row()] = float(item)
				except:
					if str(item)=="-":
						self._population_means[index.row()] = "-"
						self._population_std[index.row()] = "-"
						self._population_minmean[index.row()] = "-"
						self._population_maxmean[index.row()] = "-"
						self._population_minsd[index.row()] = "-"
						self._population_maxsd[index.row()] = "-"

				self._update_populations()


		# standard deviation of division time
		elif index.column()==3:
			#item, ok  = QtGui.QInputDialog.getDouble(self, "Change division time", "Enter standard deviation of division time:", float(current_value), 0)
			item, ok  = QtGui.QInputDialog.getText(self, "Change division time", "Enter standard deviation of division time:", QtGui.QLineEdit.Normal, current_value)
			if ok and item:
				try: 
					self._population_std[index.row()] = float(item)
				except:
					if str(item)=="-":
						self._population_means[index.row()] = "-"
						self._population_std[index.row()] = "-"
						self._population_minmean[index.row()] = "-"
						self._population_maxmean[index.row()] = "-"
						self._population_minsd[index.row()] = "-"
						self._population_maxsd[index.row()] = "-"

				self._update_populations()

		# min div time
		elif index.column()==4:
			item, ok  = QtGui.QInputDialog.getDouble(self, "Search space", "Enter minimum value for mean division time:", float(current_value), 0)
			if ok and item:
				self._population_minmean[index.row()] = item
				self._update_populations()

		# max div time
		elif index.column()==5:
			item, ok  = QtGui.QInputDialog.getDouble(self, "Search space", "Enter maximum value for mean division time:", float(current_value), 0)
			if ok and item:
				self._population_maxmean[index.row()] = item
				self._update_populations()

		# min std
		elif index.column()==6:
			item, ok  = QtGui.QInputDialog.getDouble(self, "Search space", "Enter minimum value for standard deviation:", float(current_value), 0, 1e6, 3)
			if ok and item:
				self._population_minsd[index.row()] = item
				self._update_populations()

		# max std
		elif index.column()==7:
			item, ok  = QtGui.QInputDialog.getDouble(self, "Search space", "Enter maximum value for standard deviation:", float(current_value),  0, 1e6, 3)
			if ok and item:
				self._population_maxsd[index.row()] = item
				self._update_populations()

		# info
		elif index.column()==8:
			item, ok  = QtGui.QInputDialog.getText(self, "Additional information", "Enter miscellanea information about the population:", QtGui.QLineEdit.Normal, str(current_value))
			if ok and item:
				self._population_info[index.row()] = item
				self._update_populations()


	def save_project(self):
		"""
			Two cases:		
						case 1: default filename
						 		save to filename
	
						case 2: no default filename
								invoke saveas
		"""

		if self._project_filename!=None:
			print " * Overwriting existing project"
			ret = self._save_project_to_file(self._project_filename)
		else:
			ret = self.save_project_as()
		if ret: 
			print " * Project saved correctly to %s" % self._project_filename
		else:
			print "ERROR: cannot save project"

	def save_project_as(self):
		path = QtGui.QFileDialog.getSaveFileName(self, 'Save project', ".", '*.prc')
		if path!="":
			self._save_project_to_file(path)
			self._project_filename = str(path)
			self._update_window_title()


	def _compact_populations(self):
		return self._get_population_data_frame().values.tolist()

	def _save_project_to_file(self, path):
		filehandler = open(path, 'w') 
		
		P = Project()
		
		P.initial_file = self._initial_histo_path
		P.target_file = self._target_histo_path
		P.validation_file = self._validation_histo_path
		
		P.project_name = str(self.projectname.text())
		P.populations = self._compact_populations()
		
		P.simulation_time = float(self.simulationtime.value() ) 
		P.simulation_validation_time = float(self.validationtime.value() ) 

		P.fluorescence_threshold = float(self.fluorescencethreshold.value())
		
		P.bins = int(self.bins.value())
		P.lowest_bin = float(self.lowerbin.value())
		P.highest_bin = float(self.higherbin.value())
		
		P.normalize_to_target = self.normtotarget.isChecked()
		P.asynchronous = self.async.isChecked()
		
		P.algorithm = "Fuzzy Self-Tuning PSO"
		P.swarm_size = int(self.swarmsize.value())
		P.iterations = int(self.iterations.value())		
		
		try:
			with open(path, 'wb') as output:
				pickle.dump(P, output)
				return True
		except:
			return False

	def open_project(self):
		fname = QtGui.QFileDialog.getOpenFileName(self, 'Open project', '.' ,"ProCell file (*.prc)")		
		if fname is not None:
			self._load_project_from_file(fname)
			

	def _update_window_title(self):
		self.setWindowTitle( "ProCell v" + self._version + " - " + str(self._project_filename) )

	def _reset_populations(self):
		self._population_names = []
		self._population_proportions = []
		self._population_means = []
		self._population_std = []
		self._population_minmean = []
		self._population_maxmean = []
		self._population_minsd = []
		self._population_maxsd = []
		self._population_info = []

	def _import_populations(self, pops):
		self._reset_populations()
		for p in pops:
			self._population_names.append(p[0])
			self._population_proportions.append(p[1])
			self._population_means.append(p[2])
			self._population_std.append(p[3])
			self._population_minmean.append(p[4])
			self._population_maxmean.append(p[5])
			self._population_minsd.append(p[6])
			self._population_maxsd.append(p[7])
			self._population_info.append(p[8])
		self._update_populations()

	def _load_project_from_file(self, path):
		try: 
			with open(path, 'rb') as inpt:
				P = pickle.load(inpt)
			print " * Project '%s' loaded from %s" % (P.project_name, path)
			if P is None:
				print "ERROR parsing the project file %s" % fname
				return

			self._simulated_histo = None
			self._validation_histo = None

			self._import_initial_histo(P.initial_file)
			self._update_statusbar(20)
			self._import_target_histo(P.target_file)
			self._update_statusbar(40)
			self._import_validation_histo(P.validation_file)
			self._update_statusbar(60)
			self.projectname.setText(P.project_name)
			self._import_populations(P.populations)
			self._update_statusbar(70)
			self.simulationtime.setValue(P.simulation_time)
			self.validationtime.setValue(P.simulation_validation_time)
			self.fluorescencethreshold.setValue(P.fluorescence_threshold)
			self.bins.setValue(P.bins)
			self.lowerbin.setValue(P.lowest_bin)
			self.higherbin.setValue(P.highest_bin)
			self.normtotarget.setChecked(P.normalize_to_target)
			self.async.setChecked(P.asynchronous)
			self._update_statusbar(80)
			self.swarmsize.setValue(P.swarm_size)
			self.iterations.setValue(P.iterations)
			self._update_statusbar(90)

			self.progress.reset()
			self._project_filename = str(fname)
			self._update_window_title()
		except:
			return None
		

		
	def _update_statusbar(self, v):
		self.progress.setValue(v)	

	"""
	def close(self):
		print "Shutting down ProCell"
	"""


	def _done_optimization(self):
		result_optimization = self.OPTTHREAD._dictionaries_results
		
		print " * Optimization completed"
		props, means, stdivs = result_optimization
		print props
		print means 
		print stdivs

		namepos_dict = {}
		for name, n in zip(self._population_names, range(len(self._population_names))):
			namepos_dict[name]=n

		for k,v in namepos_dict.items():
			self._population_proportions[v] = props[k]
			self._population_means[v] = means[k]
			self._population_std[v] = stdivs[k]

		self._update_populations()


	def optimize(self):
		if not self._ready_to_simulate():	
			if not self._query_for_initial():
				return
		if not self._ready_to_optimize():
			if not self._query_for_target():
				return

		try:
			os.mkdir("temp")
		except:
			print "WARNING: cannot create temporary directory (it already exists, perhaps?)"

		self._calibrator.set_types(self._population_names)
		self._calibrator.set_time_max(float(self.simulationtime.value())) #hours
		self._calibrator.set_initial_histogram_from_file(str(self._initial_histo_path))
		self._calibrator.set_target_from_file(str(self._target_histo_path))
		self._calibrator.set_output_dir("temp")
		self._calibrator.set_model_name(self.projectname.text())

		print " * All information set, creating thread"

		self.OPTTHREAD = OptimizationThread(self)
		self.connect(self.OPTTHREAD, QtCore.SIGNAL("finished()"), self._done_optimization)

		print "   Optimizer ready, launching optimization"

		self.OPTTHREAD.start()


class OptimizationThread(QThread):

	countChanged = pyqtSignal(float)

	def __init__(self, parent):
		QThread.__init__(self)
		self._parent=parent

		self._solution_optimization = None
		self._fitness_solution_optimization = None
		self._dictionaries_results = None

		self._parent.statusBar.showMessage("Optimization is starting. This process is time consuming, please wait...")
		self.timer = QtCore.QTimer(self)
		self.timer.setInterval(100)         
		self.countChanged.connect(self._parent._update_statusbar)
		self.timer.timeout.connect(self._update_status)
		self.timer.start()

	def __del__(self):
		self.timer.stop()
		self.quit()
		self.wait()

	def _update_status(self):
		pass # TODO

	def run(self):

		search_space = [ [1e-10,1.] for _ in self._parent._population_names[1:] ]

		quiescent_indices = []
		for n,v in enumerate(self._parent._population_means):
			if v!="-":
				quiescent_indices.append(n)

		copy_minmean = self._parent._population_minmean; copy_minmean = [copy_minmean[i] for i in quiescent_indices]
		copy_maxmean = self._parent._population_maxmean; copy_maxmean = [copy_maxmean[i] for i in quiescent_indices]
		copy_minsd = self._parent._population_minsd; copy_minsd = [copy_minsd[i] for i in quiescent_indices]
		copy_maxsd = self._parent._population_maxsd; copy_maxsd = [copy_maxsd[i] for i in quiescent_indices]

		search_space += [ [x,y] for (x,y) in zip(copy_minmean, copy_maxmean) ]
		search_space += [ [x,y] for (x,y) in zip(copy_minsd, copy_maxsd) ]
		#search_space += [ [x,y] for (x,y) in zip(self._population_minmean, self._population_maxmean) ]
		#search_space += [ [x,y] for (x,y) in zip(self._population_minsd, self._population_maxsd) ]


		for rep in xrange(int(self._parent.repetitions.value())):
			print " * OPTIMIZATION %d IS STARTING (please be patient...)" % (rep)
			self._solution_optimization, self._fitness_solution_optimization = \
				self._parent._calibrator.calibrate_gui(
					max_iter=int(self._parent.iterations.value()),
					swarm_size=int(self._parent.swarmsize.value()),
					search_space=search_space,
					form=self._parent,
					repetition=rep
				)

		#print self._solution_optimization, self._fitness_solution_optimization
		from estimator import fitness_gui
		self._dictionaries_results = fitness_gui(self._solution_optimization.X, arguments = {'form': self._parent}, return_dictionaries=True)
		


	def stop(self):
		print " * Trying to abort optimization..."
		self._parent._calibrator._abort_variable=True


class SimulationThread(QThread):

	countChanged = pyqtSignal(float)


	def __init__(self, parent):
		QThread.__init__(self)
		self._parent=parent
		self._parent.statusBar.showMessage("Simulation started...")
		self.result_simulation = None
		self.result_simulation_types = None
		self._what=None
		self.timer = QtCore.QTimer(self)
		self.timer.setInterval(100)         
		self.countChanged.connect(self._parent._update_statusbar)
		self.timer.timeout.connect(self._update_status)
		self._cells_in_stack = 0
		self.timer.start()

	def __del__(self):
		self.timer.stop()
		self.quit()
		self.wait()

	def get_stack_size(self):
		if self._parent.Simulator.stack is None: 
			return 0
		else:
			return self._parent.Simulator.stack.size()

	def _update_status(self):
		try:
			elapsed = 1.-1.*self.get_stack_size()/self._parent.Simulator._initial_cells_in_stack
		except:
			elapsed = 0
		self.countChanged.emit(elapsed*100)

	def _prepare_files_for_GPU(self, names, prop, mean, st):
		# each cell population = one row
		# proportion [space] mean [space] std
		with open("__temporary__", "w") as fo:
			for name in names:
				fo.write(str(prop[name])+" ")
				fo.write(str(mean[name])+" ")
				fo.write(str(st[name])+"\n")
		
	def _launch_GPU_simulation(self, executable, initial_histo, model_file, time_max, PHI, names):
		from subprocess import check_output
		from collections import defaultdict

		ret = check_output([executable, "-h", initial_histo, "-c", model_file, "-t", str(time_max), "-p", str(PHI), "-r"])
		split_rows = ret.split("\n")

		result = defaultdict(dummy)
		types  = defaultdict(list)

		for row in split_rows:
			row = row.strip("\r")
			try:
				tokenized_row = map(float, row.split("\t"))
			except ValueError:
				continue
			fluorescence = tokenized_row[0]
			total = tokenized_row[1]
			result[fluorescence]=total

			for n, amount in enumerate(tokenized_row[2:]):
				types[fluorescence].append([names[n]]*int(amount))

		return result, types


	def run(self):
		
		if self._what is None:
			print "WARNING: cannot understand what to simulate"
			return 

		if self._what=="target":
			time_max = float(self._parent.simulationtime.value() ) 
		else:
			time_max = float(self._parent.validationtime.value() ) 

		PHI      = float(self._parent.fluorescencethreshold.value())


		if self._parent._path_to_GPU_procell is not None:
			print " * Launching GPU-powered simulation"
			print "   preparing files...",

			fixed_means = map(lambda x: x if isinstance(x, float) else -1., self._parent._population_means)
			fixed_std   = map(lambda x: x if isinstance(x, float) else -1., self._parent._population_std)
			
			proportions 	= dict(zip(self._parent._population_names, self._parent._population_proportions))
			means 			= dict(zip(self._parent._population_names, fixed_means))
			stdev 			= dict(zip(self._parent._population_names, fixed_std))

			self._prepare_files_for_GPU(self._parent._population_names, proportions, means, stdev)
			self.result_simulation, self.result_simulation_types = self._launch_GPU_simulation(
				self._parent._path_to_GPU_procell, 
				self._parent._initial_histo_path, "__temporary__", 
				time_max, 
				PHI,
				self._parent._population_names
				)

			print "done"

			return
		else:

			fixed_means = map(lambda x: x if isinstance(x, float) else sys.float_info.max, self._parent._population_means)
			fixed_std   = map(lambda x: x if isinstance(x, float) else 0, self._parent._population_std)
			
			proportions 	= dict(zip(self._parent._population_names, self._parent._population_proportions))
			means 			= dict(zip(self._parent._population_names, fixed_means))
			stdev 			= dict(zip(self._parent._population_names, fixed_std))

			self.result_simulation, self.result_simulation_types = self._parent.Simulator.simulate(
				path=self._parent._initial_histo_path, 
				types=self._parent._population_names, 
				proportions=proportions,
				div_mean=means, 
				div_std=stdev, 
				time_max=time_max, 
				verbose=False, 
				phi=PHI, 
				distribution="gauss", 
				synchronous_start=False)


			if self._parent.Simulator._abort_variable:
				print " * Simulation aborted successfully"
			else:
				#print " * Simulation completed successfully"
				pass


	def stop(self):
		print " * Trying to abort simulation..."
		self._parent.Simulator._abort_variable=True


if __name__ == '__main__':
	
	app    = QtGui.QApplication(sys.argv)
	window = MainWindow()
	sys.exit(app.exec_())
