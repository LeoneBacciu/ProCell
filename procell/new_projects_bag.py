class Projects(object):

	def __init__(self):
		self._projects = []

	def add(self, p):
		if (p not in self._projects) and (p != ""):
			self._projects.append(p)
		return self._projects