from numpy import loadtxt, array
from random import random
from numpy.random import normal, gamma
from cell import Cell
import sys

def truncated_normal(mu, sigma):
	if mu==sys.float_info.max: return mu
	deviate = -1
	while(deviate<=0):
		#print mu, sigma
		deviate = normal(mu, sigma)
	return deviate


def histogram_from_list(freqlist):
	A = { key: value for key,value in freqlist }
	return A


class Stack(object):

	def __init__(self):
		self._stack = []
		self.types = []

	def push(self, C):
		self._stack.append(C)

	def pop(self):
		return self._stack.pop()

	def size(self):
		return len(self._stack)

	def _determine_type(self, types, proportions):
		if len(types)!=len(proportions):
			raise Exception("Number of types different wrt number of proportions, aborting")
		a = random()
		temp = 0
		n = 0
		while(temp<a):
			type = types[types.index(types[n])]
			#print "Checking type", type
			temp+=proportions[type]
			n+=1
		#print "Selected type", types[n-1]
		return types[n-1]

	def create_stack_from_histogram(self, 
		H0=None, types=None, proportions=None,
		div_mean=0, div_std=0, verbose=False,
		synchronous_start=True,
		distribution="gauss"
		):


		if abs(sum(proportions.values())-1)>1e-6:
			print "ERROR: proportions do not sum to 1 (%f)" % ( sum(proportions.values()) )
			exit()
		if types==None:
			print "ERROR: specify cell types"
			exit()
		if len(types)==0:
			print "ERROR: specify cell types"
			exit()
		self.types = types
		for fluo, freq in H0:
			for _ in xrange(int(freq)):
				cell = Cell(fluo)
				type = self._determine_type(types, proportions)
				cell.type = type
				cell.div_mean = div_mean[type]
				cell.div_std  = div_std[type]
				if distribution=="gauss":
					cell.timer = truncated_normal(cell.div_mean, cell.div_std)
				elif distribution=="gamma":
					if cell.div_std == 0:
						cell.timer = sys.float_info.max
					else:	
						cell.timer = gamma(cell.div_mean, cell.div_std)
				else:
					raise Exception("Distribution %s not supported" % distribution)

				# asynchronous start
				if not synchronous_start:
					if distribution=="gauss":
						cell.t = truncated_normal(cell.div_mean, cell.div_std) * random()
					elif distribution=="gamma":
						if cell.div_std==0: 
							cell.t = 0
						else:
							cell.t = gamma(cell.div_mean, cell.div_std) * random()

				self._stack.append(cell)
		if verbose: print self._stack


	def __repr__(self):
		string = ""
		for type in self.types:
			a = filter(lambda x: x.type==type, self._stack)
			string += type + ":" + str(len(a)) +" "
		return string


	def load_histogram(self, path):
		histo = loadtxt(path)
		ret = []
		n = 0
		for f,v in histo:
			if v>0:
				return histo[n:]
			n+=1
		print "ERROR: histogram empty"
		exit()


	def not_empty(self):
		return len(self._stack)>0


if __name__ == '__main__':
	
	pass