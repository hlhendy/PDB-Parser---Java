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
		if math.abs(i) > largest:
			largest = i
		elif math.abs(i) > second_largest:
			second_largest = i
	return_values.append(eig_sum*eig_sum)
	#Second largest |eigenvalue| <------------CHECK WHETHER ABS VAL
	return_values.append(second_largest)
	#Number distict eigenvalues <-------------NEED TO MAKE DISTICT VALS
	return_values.append(len(e))
	#Spectral Radius: largest |eigenvalue|
	return_values.append(largest)
	return return_values
	
	
def labelAttributes(graph):
	return_values = []
	#Link Impurity: number of edges {u, v} where L(u) != L(v) / number of edges
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
	return return_values
	
def clusterAttributes(graph):
	return_values = []
	#Closeness centrality: mean shortest path length
	return_values.append(closeness_centrality(graph))
	#Average clustering coefficient
	return_values.append(clustering(graph))
	#Small-worldness criteria
	 ##generate random graph and get closeness/clustering
	 ##compare with actual graph
	return return_values
