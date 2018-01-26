class Cell(object):

	def __init__(self, fluorescence=0, cell_type=None, timer=None, global_time=0, ID=0, verbose=False):
		if verbose:
			print " * Creating cell with fluorescence", fluorescence, "of type", cell_type, "with timer", timer
		self.fluorescence = fluorescence	
		self.type = cell_type
		self.timer = timer
		self.ID = ID
		self.T = 0.0
		self.global_time = 0.0

	def __repr__(self):
		return "CELL "+str(self.ID)+" INFO - time: "+str(self.T)+", fluorescence: " +str(self.fluorescence)

	
if __name__ == '__main__':
	
	C = Cell(fluorescence=10, timer=10)
