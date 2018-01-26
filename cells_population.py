from cell import *
from numpy import *

class CellsPopulation(object):

	def __init__(
		self,
		initial_histogram = None,
		types_p={'slow': .9, 'fast': 0.08, 'quiet': 0.02},
		timers={'slow': 1, 'fast':0.1, 'quiet': -1},
		sigma={'slow': 1, 'fast':1},
		verbose=False
	):

		self.total_cells = 0

		if abs(sum(types_p.values())-1.0)>1e-15:
			print "ERROR: probability distribution of cell types does not sum to 1, aborting."
			#print sum(types_p.values())-1.0
			exit(1)
		else:
			self.types_p = types_p

		if timers['slow']<timers['fast']:
			print "WARNING: timer of fast reactions is smaller than slow reactions."
			exit(4)
		else:
			self.timer = timers
			self.sigma = sigma

		if initial_histogram==None:
			print "ERROR: no initial histogram was provided, aborting."
			exit(2)
		else:
			self.fluorescences, self.frequencies = self.load_histogram_from_file(initial_histogram)

		self.last_ID = self.create_cells_stack()

		if verbose: print " * New population of cells created with the following characteristics:\n", self.timer, self.sigma, self.types_p


	def create_cell(self, cell_type=None, fluorescence=0, ID=0):
		if cell_type==None:
			rnd = random.random()
			if   rnd<self.types_p['slow']:
				cell_type='slow'
			elif rnd<self.types_p['slow']+self.types_p['fast']:
				cell_type='fast'
			else:
				cell_type='quiet'
		timer = self.timer[cell_type]
		if cell_type!='quiet':
			while (1):
				putative = random.normal(self.timer[cell_type], self.sigma[cell_type])	 				
				if putative>0: break
			timer = putative
		nc = Cell(fluorescence=fluorescence, timer=timer, cell_type=cell_type, ID=ID)
		return nc


	def create_cells_stack(self):
		self.stack = []
		total = 0
		for n,fluorescence in enumerate(self.fluorescences):
			for m in xrange(self.frequencies[n]):
				total += 1
				nc = self.create_cell(cell_type=None, fluorescence=fluorescence, ID=str(total))
				self.stack.append(nc)
		return total


	def load_histogram_from_file(self, path):
		A = loadtxt(path, skiprows=7)
		self.total_cells = int(sum(A.T[1]))
		return A.T[0], map(int, A.T[1])


	def filter_below_threshold(self,threshold=0):
		pos=0
		cut=0
		while(pos<len(self.fluorescences)):
			if self.fluorescences[pos]<threshold:
				pos+=1
				cut=pos
			else:
				break
		self.fluorescences=self.fluorescences[cut:]
		self.frequencies  =self.frequencies[cut:]
