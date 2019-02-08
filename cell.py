class Cell(object):

	def __init__(self, fluorescence):
		self.t = 0
		self.timer = 0
		self.type = None
		self.div_mean = 0
		self.div_std  = 0
		self.fluorescence = fluorescence

	def __repr__(self):
		return "<Cell of type=%s flu=%f>" % (self.type, self.fluorescence) 
