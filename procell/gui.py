#from __future__ import print_function
#import sip
#sip.setapi('QVariant', 2)
import pickle, os, sys
from os import sep
from PyQt5 import QtGui, QtCore, uic, QtWidgets
from numpy import loadtxt
from numpy.random import uniform
from .procell_core import Simulator
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pylab import *
import pandas as pd
#from PyQt4.QtCore import QThread, pyqtSignal
from PyQt5.QtCore import QThread, pyqtSignal
try:
	from . import resources
except:
	from . import resources3
from .estimator import hellinger1
from .project import Project
from .estimator import Calibrator, fitness_gui

from collections import defaultdict
#from PyQt4.QtGui import QApplication
from PyQt5.QtWidgets import QApplication, QMainWindow, QDialog
try:
	import seaborn as sns  
	sns.despine()
except:
	print ("WARNING: Seaborn not installed")
try:
	from configparser import ConfigParser
except:
	from ConfigParser import ConfigParser
from .new_projects_bag import Projects

try:
	MAXINT = sys.maxint
except:
	MAXINT = sys.maxsize

try:
	urange = xrange
except:
	urange = range

__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))) + sep

def dummy(): return 0

matplotlib.rcParams.update({'font.size': 8})

def spawn():
	from subprocess import check_output
	new_exec = [sys.executable, os.path.realpath(__file__), "empty"]
	ret = check_output(new_exec)

def rebin(series, lower, upper, thr=0, N=1000):
	bins = logspace(ceil(log10(lower)), ceil(log10(upper)), N)
	res  = zeros(N, dtype=int)
	pos  = 0
	for k,v in zip(series.T[0], series.T[1]):
		if k<lower or k<thr:
			if v>0: continue
		while(k>bins[pos]):
			pos+=1
			if pos>=N: 
				print ("WARNING: some fluorescence levels were above the upper limit (%f)." % k)
				return res, bins
		res[pos] += v
	return res, bins

def resampling(series, events):
	from random import sample
	total_simulated_events = int(sum(series.T[1]))
	if events>total_simulated_events:
		events_to_keep = range(total_simulated_events)
	else:
		events_to_keep = sample(range(total_simulated_events), events)
	new_dict = defaultdict(dummy)
	for etk in events_to_keep:
		#print ("_determining fluorescence for event number %d" % etk)
		accum = 0
		for k,v in zip(series.T[0], series.T[1]):
			accum += v
			if accum>=etk:
				new_dict[k] += 1
				break

	dictlist = []
	for key, value in new_dict.items():
		temp = [key,value]
		dictlist.append(temp)
	dictlist = sorted(dictlist)
	#print (dictlist)
	dictlist = array(dictlist)
	return dictlist

def verticalResizeTableViewToContents(tableView):
	count=tableView.verticalHeader().count()
	scrollBarHeight=tableView.horizontalScrollBar().height()
	horizontalHeaderHeight=tableView.horizontalHeader().height()
	rowTotalHeight = 0
	for i in urange(count):
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
			elif role == QtCore.Qt.TextAlignmentRole:
				return QtCore.Qt.AlignCenter

			if role == QtCore.Qt.BackgroundRole and index.column() == 1:
				
				if (sum(self._data["Proportion"])-1)<1e-3:
					return QtGui.QColor("#1B212F")
				else: 	
					return QtGui.QColor("#9B212F")

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
		return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable #| QtCore.Qt.ItemIsEditable

class AboutWindow(QMainWindow):

	def __init__(self):
		super(AboutWindow, self).__init__()
		uic.loadUi(__location__ +'about2.ui', self)
		
class PreferencesWindow(QDialog):

	def __init__(self):
		super(PreferencesWindow, self).__init__()
		uic.loadUi(__location__ +'preferences.ui', self)

	def browse(self):
		fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open executable file', '.', "Executable files (*.*)")
		self.path_cuprocell.setText(str(fname))

class MainWindow(QMainWindow):

	def __init__(self):
		super(MainWindow, self).__init__()
		uic.loadUi(__location__ + 'mainwindow2.ui', self)

		self._about_window = AboutWindow()
		self._preferences_window = PreferencesWindow()

		self._initial_histo_path = None
		self._target_histo_path = None
		self._validation_histo_path = None

		self._initial_histo_figure	= plt.figure(figsize=(10,7), tight_layout=True)
		self._target_histo_figure	 = plt.figure(figsize=(10,7), tight_layout=True)
		self._validation_histo_figure = plt.figure(figsize=(10,7), tight_layout=True)

		self._initial_histo_ax	= self._initial_histo_figure.add_subplot(111)
		self._target_histo_ax	 = self._target_histo_figure.add_subplot(111)
		self._validation_histo_ax = self._validation_histo_figure.add_subplot(111)

		self._initial_histo_canvas	  = FigureCanvas(self._initial_histo_figure)
		self._target_histo_canvas	  = FigureCanvas(self._target_histo_figure)
		self._validation_histo_canvas = FigureCanvas(self._validation_histo_figure)
	
		self._initial_histo_canvas.setStyleSheet("background:transparent;"); self._initial_histo_figure.set_facecolor("#ffff0000")
		self._target_histo_canvas.setStyleSheet("background:transparent;"); self._target_histo_figure.set_facecolor("#ffff0000")
		self._validation_histo_canvas.setStyleSheet("background:transparent;"); self._validation_histo_figure.set_facecolor("#ffff0000")

		"""
		self._target_histo_ax.set_xlabel("Fluorescence")
		self._target_histo_ax.set_ylabel("Cells frequency")

		self._initial_histo_ax.set_xlabel("Fluorescence")
		self._initial_histo_ax.set_ylabel("Cells frequency")

		self._validation_histo_ax.set_xlabel("Fluorescence")
		self._validation_histo_ax.set_ylabel("Cells frequency")
		"""


		self.initial_layout.addWidget(self._initial_histo_canvas)
		self.target_layout.addWidget(self._target_histo_canvas)
		self.validation_layout.addWidget(self._validation_histo_canvas)

		self._wipe_populations()
		self._wipe_histograms()

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
		self.statusBar = QtWidgets.QStatusBar()
		self.setStatusBar(self.statusBar)
		self.statusBar.showMessage("ProCell ready.")

		# progress bar
		self.progress = QtWidgets.QProgressBar()
		self.progress.setRange(0,100)
		self.statusBar.addPermanentWidget(self.progress)

		# recent menu
		self._recentmenus = self.menu_File.addMenu("&Reopen recent project")

		"""
		# project stuff
		#self._load_project_from_file("./prova.prc")		
		"""
		self._project_filename = None
		import pkg_resources
		vrs = pkg_resources.get_distribution('procell').version 

		self._version = str(vrs)

		# functionalities override
		self._force_resample = False
		self._force_CPU      = False
		self._use_dark_skin  = True

		self._unsaved_changes = False

		self._update_window_title()

		self._calibrator = Calibrator()
		self._calibrator.distribution = "gauss"  # default semantics

		self._recent_projects = Projects()
		self._config = ConfigParser()
		self._path_to_GPU_procell = None

		self._last_figure_export_path = None
		self._last_histogram_import_path = None 

		self._open_config()
		
		#self._load_project_from_file("../models/model3.prc")
		self._reload_qss()
		
		self.populations_table.verticalHeader().setVisible(False)

		self._place_logo()
		
		self.show()


	def _place_logo(self):
		label = QtWidgets.QLabel("Prova", self)
		mypix = QtGui.QPixmap (":/buttons/procell_logo_small.png")
		label.setPixmap(mypix)
		label.setSizePolicy(QtWidgets.QSizePolicy.Ignored, QtWidgets.QSizePolicy.Ignored)
		self.menubar.setCornerWidget(label, QtCore.Qt.TopRightCorner)


	def _reload_qss(self):
		sshFile=__location__ +"procell.qss"
		with open(sshFile,"r") as fh:
			self.setStyleSheet(fh.read())
		# print (" * Style Sheet %s loaded" % sshFile)

	def closeEvent(self, event):
		self._save_config()


	def _open_config(self):
		print (" * Opening last configuration")

		last_project = None

		if not os.path.exists("config.ini"):
			self._save_config()
		else:
			try:
				self._config.read('config.ini')
			except:
				print ("WARNING: config.ini not found, rebuilding...")
				self._save_config()
				return

			try: 
				ptp = self._config.get("main", "path_cuprocell")
				if ptp.strip()=="": ptp = None
				self._path_to_GPU_procell = ptp
			except:
				print ("WARNING: path to cuProCell not found")

			try:
				last_project = self._config.get("main", "last_project")
			except:
				print ("WARNING: path to last opened project not found")
			
			try:
				self._last_histogram_import_path  = self._config.get("main", "last_histo_path")
			except:
				print ("WARNING: path to last folder with histograms not found")

			try:
				self._last_figure_export_path = self._config.get("main", "last_figure_path")
			except:
				print ("WARNING: cannot determine in which directory the recent figures were exported")

			try:
				for x in urange(10):
					ppath = self._config.get("recent", "proj%d" % (x+1))
					self._recent_projects.add(ppath)
			except:
				print ("WARNING: could not retrieve recent projects list")

		self._populate_last_projects()
		self._load_project_from_file(last_project)


	def _load_old_project(self, n):
		self._load_project_from_file(self._recent_projects._projects[n])


	def _populate_last_projects(self):
		self._recentmenus.clear()
		actions = []
		for n,x in enumerate(self._recent_projects._projects):
			if x!="":
				actions.append( self._recentmenus.addAction(x) )
				if n==0: actions[-1].triggered.connect(lambda: self._load_old_project(0))
				if n==1: actions[-1].triggered.connect(lambda: self._load_old_project(1))
				if n==2: actions[-1].triggered.connect(lambda: self._load_old_project(2))
				if n==3: actions[-1].triggered.connect(lambda: self._load_old_project(3))
				if n==4: actions[-1].triggered.connect(lambda: self._load_old_project(4))
				if n==5: actions[-1].triggered.connect(lambda: self._load_old_project(5))
				if n==6: actions[-1].triggered.connect(lambda: self._load_old_project(6))
				if n==7: actions[-1].triggered.connect(lambda: self._load_old_project(7))
				if n==8: actions[-1].triggered.connect(lambda: self._load_old_project(8))
				if n==9: actions[-1].triggered.connect(lambda: self._load_old_project(9))



	def _save_config(self):
		self._config = ConfigParser()
		
		# main section
		self._config.add_section('main')
		
		if self._project_filename is None:
			self._config.set('main', 'last_project', "")
		else:
			self._config.set('main', 'last_project', str(self._project_filename))
		
		if self._path_to_GPU_procell is None:
			self._config.set('main', 'path_cuprocell', "")
		else:
			self._config.set('main', 'path_cuprocell', str(self._path_to_GPU_procell))
		
		if self._last_figure_export_path is None:
			self._config.set('main', 'last_figure_path', "")
		else:
			self._config.set('main', 'last_figure_path', str(self._last_figure_export_path))

		if self._last_histogram_import_path is None:
			self._config.set('main', 'last_histo_path', "")
		else:
			self._config.set('main', 'last_histo_path', str(self._last_histogram_import_path))

		self._config.add_section('secondary')
		self._config.set('secondary', 'force_CPU', str(self._force_CPU))
		self._config.set('secondary', 'force_resample', str(self._force_resample))

		# recent projects
		self._config.add_section('recent')
		projects = self._recent_projects._projects[-10:]
		for n in urange(10):
			try:
				pname = projects[n]
				self._config.set('recent', 'proj%d' % (n+1), pname)
			except:
				self._config.set('recent', 'proj%d' % (n+1), "")
		
		with open('config.ini', 'w') as f:
			self._config.write(f)
	

	def _open_preferences(self):
		self._preferences_window.path_cuprocell.setText(self._path_to_GPU_procell)
		self._preferences_window.forcecpu.setChecked(self._force_CPU)
		self._preferences_window.forceresampling.setChecked(self._force_resample)

		if self._preferences_window.exec_() == QtWidgets.QDialog.Accepted:
			self._path_to_GPU_procell = str(self._preferences_window.path_cuprocell.text())
			self._force_CPU = self._preferences_window.forcecpu.isChecked()
			self._force_resample = self._preferences_window.forceresampling.isChecked()

		print (" * New path to cuProCell:", self._path_to_GPU_procell)


	def _new_project(self):
		spawn()


	def _wipe_populations(self):
		self._population_names 			= []
		self._population_proportions 	= []
		self._population_means 			= []
		self._population_std   			= []
		self._population_info  			= []
		self._population_minmean 		= []
		self._population_maxmean 		= []
		self._population_minsd 			= []
		self._population_maxsd 			= []

	def _wipe_histograms(self):
		self._initial_histo	= None
		self._target_histo	 = None
		self._validation_histo = None
		self._simulated_histo  = None
		self._simulated_validation_histo = None

		#self._initial_histo_canvas	= FigureCanvas(self._initial_histo_figure)
		#self.initial_layout.addWidget (self._initial_histo_canvas)
		#self._target_histo_canvas	 = FigureCanvas(self._target_histo_figure)
		#self._validation_histo_canvas = FigureCanvas(self._validation_histo_figure)


		self._update_all_plots()

	def _new_blank(self):
		self._project_filename = None
		self.projectname.clear()

		self._wipe_populations()
		self._update_populations()

		self._wipe_histograms()

		print (" * Project wiped out")



	def _show_about(self):
		self._about_window.show()

	def _update_populations(self):
		data = self._get_population_data_frame()
		self.model = PandasModel(data=data)	
		self.populations_table.setModel(self.model)
		if not self._check_proportions():
			print (" * Error with proportions detected!")
			self._set_column_proportions_color("#FF212F")
		else:
			self._set_column_proportions_color("#1B212F")


	def _set_column_proportions_color(self, rgb):
		pass
		#for i in urange(len(self._population_proportions)):
			#item = self.populations_table.model().index(i,1)
			#self.model.setData(item, QtGui.QBrush(QtCore.Qt.red), QtCore.Qt.BackgroundRole)

		

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
		df =  pd.DataFrame({	'Population name': self._population_names, 
								'Proportion': self._population_proportions, 
								'Mean division time': self._population_means, 
								'Standard deviation': self._population_std, 
								'Mean division time (minimum)': self._population_minmean,
								'Mean division time (maximum)': self._population_maxmean,
								'Standard deviation division time (minimum)': self._population_minsd,
								'Standard deviation division time (maximum)': self._population_maxsd,
								'Info': self._population_info}
						)
		df = df.reindex([
			'Population name', 
			'Proportion',
			'Mean division time', 
			'Standard deviation',
			'Mean division time (minimum)',
			'Mean division time (maximum)',
			'Standard deviation division time (minimum)',
			'Standard deviation division time (maximum)',
			'Info'
			], axis=1)
		return df

	def drop_histogram(self):
		print ("NOT IMPLEMENTED YET")

	def _message_error(self, text, title="Import error", additional=None):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Critical)
		msg.setText(text)
		msg.setWindowTitle(title)
		if additional!=None:
			msg.setInformativeText("This is additional information")
			msg.setDetailedText(additional)
		msg.setStandardButtons(QtGui.QMessageBox.Ok)
		msg.exec_()

	def _import_initial_histo(self, path):
		path = str(path)
		try:
			A = loadtxt(path)
			print (" * Attempting the import of the initial histogram %s" % path)
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
			print (" * Attempting the import of the target histogram %s" % path)
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
			print (" * Attempting the import of the validation histogram %s" % path)
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
		fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open histogram file', self._last_histogram_import_path, "Histogram files (*.fcs *.txt)")

		# import file
		if fname!="":
			res = self._import_target_histo(fname)
			self._last_histogram_import_path = fname
			self._mark_unsaved_change()



	def import_validation_histo(self):
		# open dialog
		fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open histogram file', self._last_histogram_import_path, "Histogram files (*.fcs *.txt)")

		# import file
		if fname!="":
			res = self._import_validation_histo(fname)
			self._last_histogram_import_path = fname
			self._mark_unsaved_change()


	def _update_all_plots(self):
		if self._initial_histo 		is not None: self._update_initial_plot()
		if self._target_histo  		is not None: self._update_target_plot()
		if self._validation_histo 	is not None: self._update_validation_plot()


	def click_normalize(self):
		self._update_all_plots()


	def heuristic(self):
		if self._force_resample:	return MAXINT
		else:								return 10000


	def _update_target_plot(self, skip=0):
		self._target_histo_ax.cla()
		if self._target_histo is None: return

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())
		thr = float(self.fluorescencethreshold.value())

		if self.usecolors.isChecked():
			selected_color = "cyan"
		else:
			selected_color = "#DDDDDD"

		# plot target histogram
		if self._target_histo is not None:
			res, bins = rebin(self._target_histo, lower, higher, thr, N=calcbins)
			self._target_histo_ax.bar(bins[skip:-1], res[skip:-1], width=diff(bins[skip:]),  color=selected_color, label="Target histogram", alpha=0.5, ec="black", linewidth=0.4, align="edge")

		if self.usecolors.isChecked():
			selected_color = "green"
		else:
			selected_color = "#888888"

		if self._simulated_histo is not None:
			ratio = 1.0
			if self.normtotarget.isChecked():
				if sum(res)>self.heuristic(): 
					#print (" * Using approximate resampling")
					res2, bins2 = rebin(self._simulated_histo, lower, higher, thr, N=calcbins)
					ratio = 1.*sum(res2)/sum(res)
					
				else:
					#print (" * Using exact resampling")
					resampled = resampling(self._simulated_histo, sum(res))
					res2, bins2 = rebin(resampled, lower, higher, thr, N=calcbins)
			else:
				res2, bins2 = rebin(self._simulated_histo, lower, higher, thr, N=calcbins)
			
			self.hellingertarget.setText( "%.3f"  % hellinger1(res2/ratio, res) )
			self._target_histo_ax.bar(bins2[skip:-1], res2[skip:-1]/ratio, width=diff(bins[skip:]),  color=selected_color, label="Simulation", alpha=0.5, ec="black", linewidth=0.4, align="edge")

		self._target_histo_ax.set_xscale("symlog")
		self._target_histo_ax.set_xlim(bins[0],bins[-1])
		
		if self._use_dark_skin:
			self._target_histo_ax.set_xlabel("Fluorescence", color="#A0A7B760")
			self._target_histo_ax.set_ylabel("Cells frequency", color="#A0A7B760")
			self._target_histo_ax.tick_params(axis='both', colors='#A0A7B7', labelcolor='#A0A7B7', grid_color='#A0A7B7')
		else:
			self._target_histo_ax.set_xlabel("Fluorescence", color="black")
			self._target_histo_ax.set_ylabel("Cells frequency", color="black")
			self._target_histo_ax.tick_params(axis='both', colors='black', labelcolor='black', grid_color='black')

		if self._use_dark_skin:
			self._target_histo_ax.set_facecolor('#131721')
			self._target_histo_ax.legend().set_visible(False)
		else:
			self._target_histo_ax.set_facecolor('#FFFFFF')
			self._target_histo_ax.legend()

		self._target_histo_canvas.draw()
		self._target_histo_canvas.setToolTip(self._target_histo_path)  



	def _update_validation_plot(self, skip=0):
		self._validation_histo_ax.cla()
		if self._validation_histo is None: return

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())
		thr = float(self.fluorescencethreshold.value())

		if self.usecolors.isChecked():
			first_color = "purple"
		else:
			first_color = "#DDDDDD"
				
		if self._validation_histo is not None:
			res, bins = rebin(self._validation_histo, lower, higher, thr, N=calcbins)
			self._validation_histo_ax.bar(bins[skip:-1], res[skip:-1], width=diff(bins[skip:]),  color=first_color, label="Validation histogram", alpha=0.5, linewidth=0.4, ec="black", align="edge")

		if self.usecolors.isChecked():
			selected_color = "green"
		else:
			selected_color = "#888888"

		if self._simulated_validation_histo is not None:
			ratio = 1.0
			if self.normtotarget.isChecked():
				if sum(res)>self.heuristic(): 
					#print (" * Using approximate resampling")
					res2, bins2 = rebin(self._simulated_validation_histo, lower, higher, thr, N=calcbins)
					ratio = 1.*sum(res2)/sum(res)
					
				else:
					#print (" * Using exact resampling")
					resampled = resampling(self._simulated_validation_histo, sum(res))
					res2, bins2 = rebin(resampled, lower, higher, N=calcbins)
			else:
				res2, bins2 = rebin(self._simulated_validation_histo, lower, higher, thr, N=calcbins)

			self.hellingervalidation.setText( "%.3f"  % hellinger1(res2/ratio, res) )
			self._validation_histo_ax.bar(bins2[skip:-1], res2[skip:-1]/ratio, width=diff(bins[skip:]),  color=selected_color, label="Simulation", alpha=0.5, ec="black", linewidth=0.4, align="edge")

		self._validation_histo_ax.set_xscale("symlog")
		self._validation_histo_ax.set_xlim(bins[0],bins[-1])
		#self._validation_histo_figure.legend(framealpha=1.)
		
		if self._use_dark_skin:
			self._validation_histo_ax.set_xlabel("Fluorescence", color="#A0A7B760")
			self._validation_histo_ax.set_ylabel("Cells frequency", color="#A0A7B760")
			self._validation_histo_ax.tick_params(axis='both', colors='#A0A7B7', labelcolor='#A0A7B7', grid_color='#A0A7B7')
		else:
			self._validation_histo_ax.set_xlabel("Fluorescence", color="black")
			self._validation_histo_ax.set_ylabel("Cells frequency", color="black")
			self._validation_histo_ax.tick_params(axis='both', colors='black', labelcolor='black', grid_color='black')

		#	from matplotlib.lines import Line2D
		# custom_lines = [	Line2D([0], [0], color=first_color, lw=0, marker="o" ),
		#					Line2D([0], [0], color=selected_color, lw=0, marker="o" )]
		
		# leg = self._validation_histo_figure.legend(custom_lines, ['Validation histogram', 'Simulation'], ncol=2,  framealpha=0., loc='best')
		#for text in leg.get_texts():
		#	text.set_color("#A0A7B760")

		if self._use_dark_skin:
			self._validation_histo_ax.set_facecolor('#131721')
			leg = self._validation_histo_ax.legend()
			leg.remove()
		else:
			self._validation_histo_ax.set_facecolor('#FFFFFF')
			self._validation_histo_ax.legend()

		#self._validation_histo_figure.tight_layout()

		self._validation_histo_canvas.draw()
		self._validation_histo_canvas.setToolTip(self._validation_histo_path)  


	def new_population(self):
		self._add_population("unnamed population", proportion=0, mean=0,  st=0, minimum_mean=0, maximum_mean=0, minimum_sd=0, maximum_sd=0, info="")
		self._update_populations()	


	def remove_population(self):
		index = self.populations_table.selectionModel().selectedRows()[0]
		row = index.row()
		print (" * Removing row %d..." % row)
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


	def import_initial_histo(self):
		# open dialog
		fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open histogram file', self._last_histogram_import_path,"Histogram files (*.fcs *.txt)")

		# import file
		if fname!="":
			self._import_initial_histo(fname)
			self._last_histogram_import_path = fname
			self._mark_unsaved_change()
			
	def _update_initial_plot(self):
		self._initial_histo_ax.cla()		
		if self._initial_histo is None: 
			self._initial_histo_canvas.draw()
			return

		lower = float(self.lowerbin.value())
		higher = float(self.higherbin.value())
		calcbins = int(self.bins.value())
		
		res, bins = rebin(self._initial_histo, lower, higher, N=calcbins)

		if self._use_dark_skin:
			self._initial_histo_ax.set_xlabel("Fluorescence", color="#A0A7B760")
			self._initial_histo_ax.set_ylabel("Cells frequency", color="#A0A7B760")
			self._initial_histo_ax.tick_params(axis='both', colors='#A0A7B7', labelcolor='#A0A7B7', grid_color='#A0A7B7')
		else:
			self._initial_histo_ax.set_xlabel("Fluorescence", color="black")
			self._initial_histo_ax.set_ylabel("Cells frequency", color="black")
			self._initial_histo_ax.tick_params(axis='both', colors='black', labelcolor='black', grid_color='black')

		if self.usecolors.isChecked():
			selected_color = "red"
		else:
			selected_color = "#dddddd"
		self._initial_histo_ax.bar(bins[1:-1], res[1:-1], width=diff(bins[1:]),  
				color=selected_color, label="Initial histogram", alpha=0.5, linewidth=0.4, ec="black", align="edge")
		
		self._initial_histo_ax.set_xscale("symlog")
		self._initial_histo_ax.set_xlim(bins[0],bins[-1])
	
		custom_lines = [	Line2D([0], [0], color=selected_color, lw=0, marker="o" )]
		box = self._initial_histo_ax.get_position()
		self._initial_histo_ax.set_position([box.x0, box.y0, box.width , box.height * 0.8])
		
		# leg = self._initial_histo_figure.legend(custom_lines, ['Initial histogram'], ncol=2, bbox_to_anchor=(0.4, 0.1), framealpha=0., )
		#for text in leg.get_texts():
		#	text.set_color("#A0A7B7")
	
		if self._use_dark_skin:
			self._initial_histo_ax.set_facecolor('#131721')
			self._initial_histo_ax.legend().set_visible(False)
		else:
			self._initial_histo_ax.set_facecolor('#FFFFFF')
			self._initial_histo_ax.legend()

		#sns.despine(ax=self._initial_histo_ax, left=True, bottom=True)

		self._initial_histo_canvas.draw()
		self._initial_histo_canvas.setToolTip(self._initial_histo_path)  



	def _create_new_population_interface(self):
		return QtGui.QWidget()


	def save_figures(self):
		file = str(QtWidgets.QFileDialog.getExistingDirectory(self, "Select Directory", self._last_figure_export_path))
		if file is not None:
			self._use_dark_skin = False
			self._update_all_plots()
			self._initial_histo_figure.savefig(file+os.sep+"initial_histogram.pdf", figsize=(10,7), dpi=300)
			self._target_histo_figure.savefig(file+os.sep+"target_histogram.pdf", figsize=(10,7), dpi=300)
			self._validation_histo_figure.savefig(file+os.sep+"validation_histogram.pdf", figsize=(10,7), dpi=300)
			self._use_dark_skin = True
			self._update_all_plots()
			self._last_figure_export_path = file


	def _validate_model(self):
		pass


	def _ready_to_simulate(self):
		if self._initial_histo is not None:
			if len(self._population_names)>0:
				return True
			else:
				print ("WARNING: no populations specified, aborting.")
		else:
			print ("WARNING: initial histogram not loaded, aborting.")
		return False

	def _ready_to_optimize(self):
		if self._target_histo is not None:
			if len(self._population_names)>0:
				return True
			else:
				print ("WARNING: no populations specified, aborting.")
		else:
			print ("WARNING: target histogram not loaded, aborting.")
		return False

	def _check_proportions(self):
		return  (abs(sum(self._population_proportions)-1.0))<1e-3


	def run_simulation(self):
		if self._ready_to_simulate():
			if self._check_proportions():

				self.YTN = SimulationThread(self)
				self.YTN._what="target"
				
				self.connect(self.YTN, QtCore.SIGNAL("finished()"), self._done_simulation)

				self.launch_simulation.setEnabled(False)
				self.runvalidation.setEnabled(False)
				#self.abort.setEnabled(True)
				self.YTN.start()
			else:
				self._error_proportions()
		else:
			if self._initial_histo is None:
				self._query_for_initial()
			else:
				self._error_populations()

	def _error_proportions(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Critical)
		msg.setText("The sum of proportions is not 1. Please check the proportion of cells.")
		msg.setWindowTitle("Unable to run simulation")
		ret = msg.exec_()

	def _error_populations(self):
		msg = QtGui.QMessageBox()
		msg.setIcon(QtGui.QMessageBox.Critical)
		msg.setText("Populations were not specified. Cell populations are necessary to run ProCell.")
		msg.setWindowTitle("Unable to run simulation")
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
			#self.abort.setEnabled(True)
			self.YTN.start()

	def _done_simulation(self):
		result_simulation = self.YTN.result_simulation

		if not self.Simulator._abort_variable:

			sorted_res = array(sorted([ [a,b] for (a,b) in zip(result_simulation.keys(), result_simulation.values())]))

			if self.YTN._what == "target":
				self._simulated_histo = sorted_res
				if self.normtotarget.isChecked():
					self.statusBar.showMessage("Simulation completed, normalization in progress (please wait)...")
				self._update_target_plot()
				self.statusBar.showMessage("Simulation completed.")
			else:
				self._simulated_validation_histo = sorted_res
				if self.normtotarget.isChecked():
					self.statusBar.showMessage("Simulation for validation completed, normalization in progress (please wait)...")
				self._update_validation_plot()
				self.statusBar.showMessage("Simulation for validation completed.")
		else:
			self.statusBar.showMessage("Simulation was aborted.")

		self.launch_simulation.setEnabled(True)
		self.runvalidation.setEnabled(True)
		#self.abort.setEnabled(False)
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
			item, ok  = QtWidgets.QInputDialog.getText(self, "Change population name", "Enter new name:", QtWidgets.QLineEdit.Normal, current_value )
			if ok and item:
				self._population_names[index.row()] = item
				self._update_populations()

		# proportion
		elif index.column()==1:
			localized_input_dialog = QtWidgets.QInputDialog()
			localized_input_dialog.setInputMode(QtWidgets.QInputDialog.DoubleInput)
			localized_input_dialog.setLocale(QtCore.QLocale(QtCore.QLocale.English, QtCore.QLocale.UnitedStates))
			localized_input_dialog.setDoubleValue(float(current_value))
			localized_input_dialog.setDoubleMinimum(0)
			localized_input_dialog.setDoubleMaximum(1)
			localized_input_dialog.setDoubleDecimals(3)
			localized_input_dialog.setWindowTitle("Change proportion of population")
			localized_input_dialog.setLabelText("Enter new proportion:")
			if localized_input_dialog.exec_()==QtWidgets.QDialog.Accepted:
				item = localized_input_dialog.doubleValue()
				self._population_proportions[index.row()] = item
				self._update_populations()

		# mean division time
		elif index.column()==2:
			# item, ok  = QtGui.QInputDialog.getDouble(self, "Change division time", "Enter new mean division time:", float(current_value), 0)
			item, ok  = QtWidgets.QInputDialog.getText(self, "Change division time", "Enter new mean division time:", QtWidgets.QLineEdit.Normal, current_value)
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
			#item, ok  = QtWidgets.QInputDialog.getDouble(self, "Change division time", "Enter standard deviation of division time:", float(current_value), 0)
			item, ok  = QtWidgets.QInputDialog.getText(self, "Change division time", "Enter standard deviation of division time:", QtWidgets.QLineEdit.Normal, current_value)
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
			item, ok  = QtWidgets.QInputDialog.getDouble(self, "Search space", "Enter minimum value for mean division time:", float(current_value), 0)
			if ok and item:
				self._population_minmean[index.row()] = item
				self._update_populations()

		# max div time
		elif index.column()==5:
			item, ok  = QtWidgets.QInputDialog.getDouble(self, "Search space", "Enter maximum value for mean division time:", float(current_value), 0)
			if ok and item:
				self._population_maxmean[index.row()] = item
				self._update_populations()

		# min std
		elif index.column()==6:
			item, ok  = QtWidgets.QInputDialog.getDouble(self, "Search space", "Enter minimum value for standard deviation:", float(current_value), 0, 1e6, 3)
			if ok and item:
				self._population_minsd[index.row()] = item
				self._update_populations()

		# max std
		elif index.column()==7:
			item, ok  = QtWidgets.QInputDialog.getDouble(self, "Search space", "Enter maximum value for standard deviation:", float(current_value),  0, 1e6, 3)
			if ok and item:
				self._population_maxsd[index.row()] = item
				self._update_populations()

		# info
		elif index.column()==8:
			item, ok  = QtWidgets.QInputDialog.getText(self, "Additional information", "Enter miscellanea information about the population:", QtWidgets.QLineEdit.Normal, str(current_value))
			if ok and item:
				self._population_info[index.row()] = item
				self._update_populations()

		self._mark_unsaved_change()


	def _check_proliferating_for_pe(self):
		if len(list(filter(lambda x: x!='-', self._population_minsd)))==0:

			print ("WARNING: no proliferating populations found, I cannot calibrate")

			msg = QtGui.QMessageBox()
			msg.setIcon(QtGui.QMessageBox.Critical)
			msg.setText("All the populations seem to be quiescent. Please add at least one proliferating population.")
			msg.setWindowTitle("Unable to run optimization")
			#msg.setStandardButtons(QtGui.QMessageBox.OK)
			ret = msg.exec_()

			return False
		else:
			return True


	def _check_boundaries_for_pe(self):
		all_values = self._population_minsd + self._population_maxsd +  self._population_minmean + self._population_maxmean
		ready = len(list(filter(lambda x: x==0, all_values)))==0

		if not ready:
			print ("WARNING: some boundaries are not set")

			msg = QtGui.QMessageBox()
			msg.setIcon(QtGui.QMessageBox.Critical)
			msg.setText("Some boundaries were not set properly. Please set all boundaries (for both mean and standard deviation) before running parameters fitting.")
			msg.setWindowTitle("Unable to run optimization")
			#msg.setStandardButtons(QtGui.QMessageBox.OK)
			ret = msg.exec_()

			return False
		else:
			return True


	def save_project(self):
		"""
			Two cases:		
						case 1: default filename
						 		save to filename
	
						case 2: no default filename
								invoke saveas
		"""

		if self._project_filename!=None:
			print (" * Overwriting existing project")
			ret = self._save_project_to_file(self._project_filename)
		else:
			ret = self.save_project_as()
		if ret: 
			print (" * Project saved correctly to %s" % self._project_filename)
			self._unmark_unsaved_change()
		else:
			print ("ERROR: cannot save project")

	def save_project_as(self):
		path = QtWidgets.QFileDialog.getSaveFileName(self, 'Save project', ".", '*.prc')
		if path!="":
			self._save_project_to_file(path)
			self._project_filename = str(path)
			self._update_window_title()


	def _compact_populations(self):
		return self._get_population_data_frame().values.tolist()


	def _gui_update(self):
		self._update_all_plots()
		self._mark_unsaved_change()
		self.my_timer.stop()

	def _simfield_updated(self):
		self.my_timer = QtCore.QTimer()
		self.my_timer.timeout.connect(self._gui_update)
		self.my_timer.start(1000) #1 min intervall
		

	def _save_project_to_file(self, path):
				
		P = Project()
		
		P.initial_file = self._initial_histo_path
		P.target_file = self._target_histo_path
		P.validation_file = self._validation_histo_path

		P.initial_histo = self._initial_histo
		P.target_histo = self._target_histo
		P.validation_histo = self._validation_histo
		
		P.project_name = str(self.projectname.text())
		P.populations = self._compact_populations()
		
		P.simulation_time = float(self.simulationtime.value() ) 
		P.simulation_validation_time = float(self.validationtime.value() ) 

		P.fluorescence_threshold = float(self.fluorescencethreshold.value())
		
		P.bins = int(self.bins.value())
		P.lowest_bin = float(self.lowerbin.value())
		P.highest_bin = float(self.higherbin.value())
		
		P.normalize_to_target = self.normtotarget.isChecked()
		P.asynchronous = self.asyncr.isChecked()
		
		P.algorithm = "Fuzzy Self-Tuning PSO"
		P.swarm_size = int(self.swarmsize.value())
		P.iterations = int(self.iterations.value())		

		try:
			filehandler = open(path, 'w') 
		except:
			self._message_error("Unable to write on file %s.\nFile might be write protected or locked by another application." % path, title="Error while saving project file")
			return False
		
		try:
			with open(path, 'wb') as output:
				pickle.dump(P, output)
				return True
		except:
			return False

	def open_project(self):
		fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open project', '.' ,"ProCell file (*.prc)")		
		if fname is not None:
			self._load_project_from_file(fname)
			

	def _update_window_title(self):
		if  self._unsaved_changes:
			edited = "*"
		else:
			edited = "" 
		title_string = "%s%s - ProCell %s" % ( str(self._project_filename), edited, self._version)
		self.setWindowTitle( title_string )


	def _mark_unsaved_change(self):
		self._unsaved_changes = True
		self._update_window_title()


	def _unmark_unsaved_change(self):
		self._unsaved_changes = False
		self._update_window_title()


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
		
		if path=="" or path is None: return None

		with open(path, 'rb') as inpt:
			P = pickle.load(inpt)
		

		try: 
			with open(path, 'rb') as inpt:
				P = pickle.load(inpt)
			print (" * Project '%s' loaded from %s" % (P.project_name, path))
			if P is None:
				print ("ERROR parsing the project file %s" % path)
				return

			self._simulated_histo = None
			self._validation_histo = None
		except:
			print ("Unable to open the project file")
			return None

		try:
			if P.initial_file is not None:
				self._import_initial_histo(P.initial_file)
			self._update_statusbar(20)
		except:
			print ("Unable to import the initial histogram")
			return None

		try:
			if P.target_file is not None:
				self._import_target_histo(P.target_file)
			self._update_statusbar(40)
		except:
			print ("Unable to import the target histogram")
			return None

		try:
			if P.validation_file is not None:
				self._import_validation_histo(P.validation_file)
			self._update_statusbar(60)
		except:
			print ("Unable to import the validation histogram")
			return None

		try:
			self.projectname.setText(P.project_name)
			self._import_populations(P.populations)
			self._update_statusbar(70)
		except:
			print ("Unable to import the populations")
			return None

		try:
			self.simulationtime.setValue(P.simulation_time)
			self.validationtime.setValue(P.simulation_validation_time)
			self.fluorescencethreshold.setValue(P.fluorescence_threshold)
			self.bins.setValue(P.bins)
			self.lowerbin.setValue(P.lowest_bin)
			self.higherbin.setValue(P.highest_bin)
		except:
			print ("Unable to set the time, or fluorescence, or bins")
			return None

		try:
			self.normtotarget.setChecked(P.normalize_to_target)
			self.asyncr.setChecked(P.asynchronous)
			self._update_statusbar(80)
			self.swarmsize.setValue(P.swarm_size)
			self.iterations.setValue(P.iterations)
			self._update_statusbar(90)
		except:
			print ("Unable to set the PE values")
			return None

		self.progress.reset()
		self._project_filename = str(path)
		self._recent_projects.add(self._project_filename)
		self._populate_last_projects()
		self._update_window_title()
			

		
	def _update_statusbar(self, v):
		self.progress.setValue(v)	

	"""
	def close(self):
		print "Shutting down ProCell"
	"""

	def _save_simtarget_file(self):
		if self._simulated_histo is not None:
			path = QtWidgets.QFileDialog.getSaveFileName(self, 'Save histogram', ".", '*.txt')
			if path is not None:
				savetxt(path, self._simulated_histo)


	def _done_optimization(self):
		result_optimization = self.OPTTHREAD._dictionaries_results
		
		print (" * Optimization completed")
		props, means, stdivs = result_optimization
		print (props)
		print (means)
		print (stdivs)

		namepos_dict = {}
		for name, n in zip(self._population_names, range(len(self._population_names))):
			namepos_dict[name]=n

		for k,v in namepos_dict.items():
			self._population_proportions[v] = props[k]
			self._population_means[v] = means[k]
			self._population_std[v] = stdivs[k]

		self._update_populations()
		self.run_simulation()


	def optimize(self):


		if not self._ready_to_simulate():	
			if not self._query_for_initial():
				return

		if not self._ready_to_optimize():
			if not self._query_for_target():
				return
		
		if not self._check_boundaries_for_pe():
			return

		if not self._check_proliferating_for_pe():
			return

		try:
			os.mkdir("temp")
		except:
			print ("WARNING: cannot create temporary directory (it already exists, perhaps?)")

		self._calibrator.set_types(self._population_names)
		self._calibrator.set_time_max(float(self.simulationtime.value())) #hours
		self._calibrator.set_initial_histogram_from_file(str(self._initial_histo_path))
		self._calibrator.set_target_from_file(str(self._target_histo_path))
		self._calibrator.set_output_dir("temp")
		self._calibrator.set_model_name(self.projectname.text())

		print (" * All information set, creating thread")

		self.OPTTHREAD = OptimizationThread(self)
		self.connect(self.OPTTHREAD, QtCore.SIGNAL("finished()"), self._done_optimization)

		print ("   Optimizer ready, launching optimization")

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


		for rep in urange(int(self._parent.repetitions.value())):
			print (" * OPTIMIZATION %d IS STARTING (please be patient...)" % (rep))
			self._solution_optimization, self._fitness_solution_optimization = \
				self._parent._calibrator.calibrate_gui(
					max_iter=int(self._parent.iterations.value()),
					swarm_size=int(self._parent.swarmsize.value()),
					search_space=search_space,
					form=self._parent,
					repetition=rep,
					loginit=self._parent.logarithmic.isChecked()
				)

		#print self._solution_optimization, self._fitness_solution_optimization
		self._dictionaries_results = fitness_gui(self._solution_optimization.X, arguments = {'form': self._parent}, return_dictionaries=True)
		


	def stop(self):
		print (" * Trying to abort optimization...")
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
		
		ret = check_output([executable, "-h", initial_histo, "-c", model_file, "-t", str(time_max), "-p", str(PHI), "-r"])
		ret = str(ret.decode('ascii')).replace('\n\n', '\n')
		with open('pd', 'w') as fo:
			fo.write(ret)
		split_rows = ret.split("\n")
		#print (split_rows)
		#exit()
		result = defaultdict(dummy)
		types  = defaultdict(list)

		for row in split_rows:
			row = row.strip("\r")
			try:
				tokenized_row = list(map(float, row.split("\t")))
			except ValueError:
				continue
			#print (tokenized_row)
			fluorescence = tokenized_row[0]
			total = tokenized_row[1]
			result[fluorescence]=total

			for n, amount in enumerate(tokenized_row[2:]):
				types[fluorescence].append([names[n]]*int(amount))

		return result, types


	def run(self):
		
		if self._what is None:
			print ("WARNING: cannot understand what to simulate")
			return 

		if self._what=="target":
			time_max = float(self._parent.simulationtime.value() ) 
		else:
			time_max = float(self._parent.validationtime.value() ) 

		PHI	  = float(self._parent.fluorescencethreshold.value())


		if (self._parent._path_to_GPU_procell is not None) and (not self._parent._force_CPU):
			print (" * Launching GPU-powered simulation")
			#print ("   preparing files...",)

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

			#print ("done")

			return
		else:
			print (" * Launching CPU-based simulation")

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
				print (" * Simulation aborted successfully")
			else:
				#print " * Simulation completed successfully"
				pass


	def stop(self):
		print (" * Trying to abort simulation...")
		self._parent.Simulator._abort_variable=True


def main():
	app	= QtWidgets.QApplication(sys.argv)
	window = MainWindow()
	sys.exit(app.exec_())


if __name__ == '__main__':
	
	app	= QtWidgets.QApplication(sys.argv)
	window = MainWindow()
	sys.exit(app.exec_())
