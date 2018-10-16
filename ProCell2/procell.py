from create_stack import Stack, truncated_normal
from cell import Cell
from numpy.random import normal
from collections import defaultdict
import sys


def dummy(): return 0




class Simulator(object):

	def __init__(self):
		self.stack = None

	def simulate(self, 
		path=None, 
		types=None, 
		proportions=None,	
		div_mean=None, 
		div_std=None, 
		time_max=0, 
		verbose=False,
		phi=None,
		synchronous_start=True
		):

		self.stack = Stack()
		histo = self.stack.load_histogram(path)
		self.stack.create_stack_from_histogram( 
			H0=histo, types=types, proportions=proportions, 
			div_mean=div_mean, div_std=div_std, verbose=verbose,
			synchronous_start=synchronous_start)

		result = defaultdict(dummy)
		types  = defaultdict(list)
		if phi==None: 	phi_1 = histo[0][0]
		else: 			phi_1 = phi

		while self.stack.not_empty()>0:
			curcell = self.stack.pop()
			if curcell.t + curcell.timer > time_max:
				result[curcell.fluorescence] += 1
				types[curcell.fluorescence].append(curcell.type)
			else:
				curcell.t += curcell.timer
				if curcell.fluorescence / 2 > phi_1:
					daughter1 = Cell(curcell.fluorescence / 2)
					daughter2 = Cell(curcell.fluorescence / 2)
					daughter1.t = curcell.t
					daughter2.t = curcell.t
					daughter1.type = curcell.type
					daughter2.type = curcell.type
					daughter1.timer = truncated_normal( div_mean[curcell.type], div_std[curcell.type] )
					daughter2.timer = truncated_normal( div_mean[curcell.type], div_std[curcell.type] )
					self.stack.push(daughter1)
					self.stack.push(daughter2)

		return result, types


if __name__ == '__main__':
	
	pass