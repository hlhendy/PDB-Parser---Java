import networkx as nx
import numpy
import math
import operator

def shortest_path_lengths(graph):
	n = nx.number_of_nodes(graph)
	spl = []
	for i in range(n):
		for j in range(n):
			if i != j:
				try:
					spl.append(nx.shortest_path_length(graph, i, j))
				except nx.NetworkXNoPath:
					spl.append(0)
	return spl

def generalAttributes(graph):
	return_values = []
	#Number of edges
	return_values.append(nx.number_of_edges(graph))	
	#Density
	return_values.append(nx.density(graph))
	#Average Degree
	num_endpoints = 0
	degrees = nx.degree(graph)
	deg_sum = 0
	for i in degrees:
		if degrees[i] == 1:
			num_endpoints += 1
		deg_sum += degrees[i]
	return_values.append(float(deg_sum)/len(degrees))
	#Percentage of endpoints(number of nodes with deg=1)	
	return_values.append(float(num_endpoints)/len(degrees))
	return return_values

def eccentricityAttributes(graph):
	return_values = []
	#Average effective eccentricity
	eccVals = []
	e = 0
	for n in graph.nodes():
		try: 
			eccVals.append(nx.eccentricity(graph, v=n))
		except nx.NetworkXError:
			eccVals.append(0)
	eccSum = 0
	center_nodes = 0
	diameter = max(eccVals)
	radius = min(eccVals)
	for i in range(len(eccVals)):
		if eccVals[i] ==  radius:
			center_nodes += 1
		eccSum += eccVals[i]
	return_values.append(eccSum / float(nx.number_of_nodes(graph)))	
	#Effective diameter
	return_values.append(diameter)
	#Effective radius
	return_values.append(radius)
	#Percentage central nodes
	return_values.append(center_nodes / float(nx.number_of_nodes(graph)))
	return return_values

def eigenvalueAttributes(graph):
	return_values = []
	#Compute eigenvalues on L as ndarray object
	#L = nx.normalized_laplacian_matrix(graph)
	#e = numpy.linalg.eigvals(L.A)
	e = nx.adjacency_spectrum(graph) #numPy array
	#energy: squared sum of eigenvalues
	eig_sum = 0
	largest = 0
	second_largest = 0
	unique_values = []
	for i in e:
		eig_sum += (i*i)
		if i < 0:
			abs_i = i * -1
		else:
			abs_i = i
		if abs_i > largest:
			largest = i
		elif abs_i > second_largest:
			second_largest = i
		if i not in unique_values:
			unique_values.append(i)
	return_values.append(eig_sum)
	#Second largest eigenvalue - http://docs.scipy.org/doc/numpy-1.10.1/reference/generated/numpy.linalg.eigvals.html
	## -------> returns values with multiplicity and not in any order; not necessarily real for real matrices
	return_values.append(second_largest)
	#Number distict eigenvalues
	return_values.append(len(unique_values))
	#Spectral Radius: largest |eigenvalue|
	return_values.append(largest)
	return return_values
	
	
def labelAttributes(graph):
	return_values = []
	num_edges = 0
	edges_sum = 0
	#Link Impurity: number of edges {u, v} where L(u) != L(v) / number of edges
	for u,v in graph.edges():
		num_edges += 1
		if(graph.node[u]['hydro'] != graph.node[v]['hydro']):
			edges_sum += 1
	return_values.append(float(edges_sum)/num_edges)
	#Neighborhood Impurity: vertex degree for neighboring nodes with different label
	deg = 0
	num_nodes = 0
	for cn in graph.nodes():
		num_nodes += 1
		neighborhood = graph.neighbors(cn)
		for nn in neighborhood:
			if graph.node[cn]['hydro'] != graph.node[nn]['hydro']:
				deg += 1
	return_values.append(deg/num_nodes)		
	return return_values
	
def clusterAttributes(graph):
	return_values = []
	#Closeness centrality: mean shortest path length
	#nx.closeness_centrality returns dict with val for each node
	#need to find average
	closenessVals = nx.closeness_centrality(graph)
	c = 0
	for v in closenessVals:
		c += closenessVals[v]
	c = c / len(closenessVals)
	return_values.append(c)
	#Average clustering coefficient
	clustering = nx.average_clustering(graph)
	return_values.append(clustering)
	return return_values

def smallWorldness(graph):
	return_values = []
	#Small-worldness criteria
	n = len(nx.nodes(graph))
	e = len(nx.edges(graph))
	#probability of edges: (number of edges in real graph)/possible edges
	p = e/float((n*(n-1)/2.0))
	#generate random graph using probability
	rand_graph = nx.fast_gnp_random_graph(n, p)
	#calculate values for real graph and random graph
	Creal = nx.transitivity(graph) #float
	Crand = nx.transitivity(rand_graph) #float
	Lreal = 0
	Lrand = 0
	real_sum = 0
	rand_sum = 0
	splReal = shortest_path_lengths(graph)
	splRand = shortest_path_lengths(rand_graph)
	for i in range(len(splReal)):
		real_sum += splReal[i]
		rand_sum += splRand[i]
	Lreal = real_sum / len(splReal)
	Lrand = rand_sum / len(splRand)		
	#compare with actual graph
	if(Lreal != 0 and Lrand !=0 and Crand !=0):
		S = (Creal)/(Crand) / (float(Lreal)/(Lrand))
	else:
		S = 0
	return_values.append(S)
	return return_values
