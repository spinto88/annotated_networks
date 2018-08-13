import ctypes 

class Edge(ctypes.Structure):

	_fields_ = [('target', ctypes.c_int), 
		('weight', ctypes.c_double)]

class Vertex(ctypes.Structure):

	_fields_ = [('id', ctypes.c_int),
		('degree', ctypes.c_int),
		('label', ctypes.POINTER(c_char)),
		('edge', ctypes.POINTER(Edge))]

class Network(ctypes.Structure):

	_fields_ = [('nvertices', ctypes.c_int),
		('directed', ctypes.c_int),
		('vertex', ctypes.POINTER(Vertex))]

	# Read GML file
	def readgml(self, file_name):
		pass

	# To put here belief propagation by steps
	def community_detection(self):
		pass






