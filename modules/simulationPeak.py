class Peak:
	def _init_(self, x, y):
		self.X = (x,y)
		self.V = (0,0)
		self.A = (0,0)
		self.neighbors = []