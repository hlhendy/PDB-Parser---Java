import networkx as nx
import numpy
import math

def eigenvalueAttributes(graph):
	return_values = []
	#Compute eigenvalues on L as ndarray object
	L = nx.normalized_laplacian_matrix(graph)
	e = numpy.linalg.eigvals(L.A)
	#energy: squared sum of eigenvalues
	eig_sum = 0
	largest = 0
	second_largest = 0
	for i in e:
		eig_sum += i
		if math.fabs(i) > largest:
			largest = i
		elif math.fabs(i) > second_largest:
			second_largest = i
	return_values.append(eig_sum*eig_sum)
	#Second largest |eigenvalue| <------------CHECK WHETHER ABS VAL
	return_values.append(second_largest)
	#Number distict eigenvalues
	unique_values = []
	for v in e:
		if v not in unique_values:
			unique_values.append(v)
	return_values.append(len(unique_values))
	#Spectral Radius: largest |eigenvalue|
	return_values.append(largest)
	return return_values
	
	
def labelAttributes(graph):
	return_values = []
	num_edges = 0
	edges_sum = 0
	#Link Impurity: number of edges {u, v} where L(u) != L(v) / number of edges
<<<<<<< HEAD
	for e in graph.edges_iter():
		num_edges += 1
		if(e[0] != e[1]):
			edges_sum += 1
	return_values.append(edges_sum/num_edges)
	#Neighborhood Impurity: vertex degree for only those with different label
	deg = 0
	num_nodes = 0
	for cn in graph.nodes_iter():
		num_nodes += 1
		neighborhood = graph.neighbors(cn)
		for nn in neighborhood:
			if cn != nn:
				deg += 1
	return_values.append(deg/num_nodes)
=======
	total_edges = len(graph.edges())
	noMatchEdges = 0
	#node list
	nodeList = graph.nodes(data=True)
	for i,j in enumerate(nodeList[:-1]):
		if j != nodeList[i+1]:
			noMatchEdges +=  1
	return_values.append(noMatchEdges/total_edges)		
	#Neighborhood Impurity: vertex degree for only those with different label
	deg = 0
	for node in graph:
		neighbors = graph.neighbors(node)
		for n,i in enumerate(neighbors[:-1]):
			#if labels match, increment degree
			if n == neighbors[i+1]:
				deg += 1
	return_values.append(deg/len(graph.nodes())
	#Label Entropy: For a graph G with k labels, 
	##E(G) = - Sum(prob of given label * log prob of given label)
>>>>>>> 21fede94458d13ec4aaa92ebe6d4ed71ef1245d4
	return return_values
	
def clusterAttributes(graph):
	return_values = []
	#Closeness centrality: mean shortest path length
	#nx.closeness_centrality returns dict with val for each node
	#need to find average
	LgVals = nx.closeness_centrality(graph)
	Lg = 0
	for v in LgVals:
		Lg += LgVals[v]
	Lg = Lg / len(LgVals)
	return_values.append(Lg)
	#Average clustering coefficient
	Cg = nx.average_clustering(graph)
	return_values.append(Cg)
	#Small-worldness criteria
	##generate random graph and get closeness/clustering
	n = len(nx.nodes(graph))
	e = len(nx.edges(graph))
	#probability of edges: number of edges in real graph/possible edges
	p = e/(n*(n-1)/2)
	rand_graph = nx.fast_gnp_random_graph(n, p)
	LrandVals = nx.closeness_centrality(rand_graph)
	Lrand = 0
	for v in LgVals:
		Lrand += LrandVals[v]
	Lrand = Lrand / len(LrandVals)
	Crand = nx.average_clustering(rand_graph)
	##compare with actual graph
	if(Lg != 0 and Lrand !=0 and Crand !=0):
		S = (Cg/Crand) / (Lg/Lrand)
	else:
		S = 0
	return_values.append(S)
	return return_values
