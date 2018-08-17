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

def graphtool2newman(gt_graph, label):

    # Create Newman graph to be used in C programs
    graph = Network()

    # Number of vertices
    graph.nvertices = gt_graph.num_vertices()

    # Is directed?
    graph.directed = int(gt_graph.is_directed())
  
    # Graph vertices array
    graph.vertex = (graph.nvertices * Vertex)()

    for i in range(graph.nvertices):

        # Vertex id
        graph.vertex[i].id = i 

        # Vertex degree
        graph.vertex[i].degree = gt_graph.vertex(i).out_degree()

        # Vertex edge array
        graph.vertex[i].edge = (gt_graph.vertex(i).out_degree() * Edge)()

        # Neighbors assigment 
        vs = gt_graph.vertex(i)
        neigh = [int(k) for k in vs.out_neighbors()]
        for j in range(gt_graph.vertex(i).out_degree()):
            graph.vertex[i].edge[j].target = neigh[j]

    # Label information
    labels = list(set(gt_graph.vertex_properties[label]))

    # Number of distinct labels 
    nmlabels = len(labels)

    # Integer value for each categorical label
    x = [labels.index(gt_graph.vertex_properties[label][i]) for i in range(gt_graph.num_vertices())]

    # Number of vertices woth each label
    nx = [list(gt_graph.vertex_properties[label]).count(l) for l in labels]

    # Return the graph with information, number of distinct labels, 
    # integer value per vertex, and number of vertices with each label
    return graph, nmlabels, x, nx

def core_function(graph, k_comm, nmlabels, x, nx, steps, random_seed):

    import os
    import subprocess

    make_process = subprocess.Popen("make", stderr=subprocess.STDOUT)
    make_process.wait()

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

    # Initialization of metadata 
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
    lparams = []
    while(step < steps):

        libc.bp(graph, k_comm, x, gmma, omega, eta, q)

        lparams.append(libc.params(graph, k_comm, x, nx, nmlabels, nrx, gmma,\
                omega, eta, q))

        step += 1

    communities = [[q[i][k] for k in range(k_comm)] \
                    for i in range(graph.nvertices)]

    mix_matrix = np.array([[gmma[k][l] for l in range(nmlabels)] \
                  for k in range(k_comm)])

    return communities, mix_matrix, lparams


def annotated_community_detection(gt_graph, k_comm = 2, label = 'label', steps = 100, random_seed = 123457):

    graph, nmlabels, x, nx = graphtool2newman(gt_graph, label)

    comm, mixm, l = core_function(graph, k_comm, nmlabels, x, nx, steps, random_seed)

    return comm, mixm, l
