import ctypes 
import numpy as np

class Edge(ctypes.Structure):

    _fields_ = [('target', ctypes.c_int), 
		('weight', ctypes.c_double)]

    def __init__(self):
        pass

class Vertex(ctypes.Structure):

    _fields_ = [('id', ctypes.c_int),
		('degree', ctypes.c_int),
		('label', ctypes.POINTER(ctypes.c_char)),
		('edge', ctypes.POINTER(Edge))]

    def __init__(self):
        pass

class Network(ctypes.Structure):

    _fields_ = [('nvertices', ctypes.c_int),
		('directed', ctypes.c_int),
		('vertex', ctypes.POINTER(Vertex))]

    def __init__(self):
        pass

def pynet2newman():

    from graph_tool import Graph

    g = Graph()
    g.load(file_name = 'sbm-meta.gml')

    graph = Network()

    graph.nvertices = g.num_vertices()
    graph.directed = 0
    graph.vertex = (graph.nvertices * Vertex)()

    for i in range(graph.nvertices):
        graph.vertex[i].id = i
        graph.vertex[i].degree = g.vertex(i).out_degree()
#        graph.vertex[i].label = 'Meta'

        graph.vertex[i].edge = (g.vertex(i).out_degree() * Edge)()
        vs = g.vertex(i)
        neigh = [int(k) for k in vs.out_neighbors()]
        for j in range(g.vertex(i).out_degree()):
            graph.vertex[i].edge[j].target = neigh[j]

    return graph

def getmetadata():

    from graph_tool import Graph

    g = Graph()
    g.load(file_name = 'sbm-meta.gml')

    labels = list(set(g.vertex_properties['label']))

    nmlabels = len(labels)

    x = [labels.index(g.vertex_properties['label'][i]) for i in range(g.num_vertices())]

    nx = [list(g.vertex_properties['label']).count(l) for l in labels]

     # To translate graph_tool, igraph, or networkx graphs to C structure defined by Mark Newman
    return nmlabels, x, nx

def core_function(graph, k_comm = 2, steps = 100, random_seed = 123467):

    import os
    libc = ctypes.CDLL(os.getcwd() + '/libc.so')

    libc.bp.argtypes = [Network,\
                       ctypes.c_int,\
                       ctypes.POINTER(ctypes.c_int),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
    libc.bp.restype = ctypes.c_int

    libc.params.argtypes = [Network,\
                       ctypes.c_int,\
                       ctypes.POINTER(ctypes.c_int),\
                       ctypes.POINTER(ctypes.c_int),\
                       ctypes.c_int,\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double)),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_double))),\
                       ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
    libc.params.restype = ctypes.c_double
                           
    np.random.seed(random_seed)

    # Inicializar metadata 
    nmlabels, x, nx = getmetadata()
    x = (graph.nvertices * ctypes.c_int)(*x)
    nx = (graph.nvertices * ctypes.c_int)(*nx)

    # q has information about community assigment for each node u
    # It is initializated at random    
    q = (graph.nvertices * ctypes.POINTER(ctypes.c_double))()
    for u in range(graph.nvertices):
        aux = np.random.dirichlet(np.ones(k_comm))
        q[u] = (k_comm * ctypes.c_double)(*aux)
  

    eta = (graph.nvertices * ctypes.POINTER(ctypes.POINTER(ctypes.c_double)))()
    for u in range(graph.nvertices):
        eta[u] = (graph.vertex[u].degree * ctypes.POINTER(ctypes.c_double))()
        for i in range(graph.vertex[u].degree):
            eta[u][i] = (k_comm * ctypes.c_double)()
            v = graph.vertex[u].edge[i].target
            for r in range(k_comm):
                eta[u][i][r] = q[v][r]

    nrx = (k_comm * ctypes.POINTER(ctypes.c_double))()
    for r in range(k_comm):
        nrx[r] = (nmlabels * ctypes.c_double)()

    gmma = (k_comm * ctypes.POINTER(ctypes.c_double))()
    for r in range(k_comm):
        aux = np.random.dirichlet(np.ones(nmlabels))
        gmma[r] = (nmlabels * ctypes.c_double)(*aux)
        
    # Twice the number of edges
    twom = np.sum([graph.vertex[u].degree for u in range(graph.nvertices)])

    omega = (k_comm * ctypes.POINTER(ctypes.c_double))()
    for r in range(k_comm):
        omega[r] = (k_comm * ctypes.c_double)()
        for s in range(k_comm):
            if r == s:
                omega[r][s] = np.random.random()/twom;
            elif r < s:
                omega[r][s] = np.random.random()/twom;
            else:
                omega[r][s] = omega[s][r]

    step = 0
    while(step < steps):

        libc.bp(graph, k_comm, x, gmma, omega, eta, q)

        libc.params(graph, k_comm, x, nx, nmlabels, nrx, gmma,\
                omega, eta, q)

        step += 1

    communities = [[q[i][k] for k in range(k_comm)] \
                    for i in range(graph.nvertices)]

    mix_matrix = np.array([[gmma[l][k] for l in range(nmlabels)] \
                  for k in range(k_comm)])

    return communities, mix_matrix

def main():

    graph = pynet2newman()
    comm, mixm = core_function(graph, k_comm = 4)

    print comm[:10]
    print mixm

if __name__ == "__main__":
    main()


