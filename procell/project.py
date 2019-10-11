class Project(object):
	
	"""docstring for Project"""
	def __init__(self):
		super(Project, self).__init__()		
		self.initial_file = None
		self.target_file = None
		self.validation_file = None
		
		self.project_name = None
		self.populations  = None
		
		self.simulation_time = None
		self.simulation_validation_time = None
		
		self.fluorescence_threshold = None
		
		self.bins = None
		self.lowest_bin = None
		self.highest_bin = None
		
		self.normalize_to_target = None
		self.asynchronous = None
		
		self.algorithm = None
		self.swarm_size = None
		self.iterations = None