def eigenvalueAttributes(graph):
	return_values = []
	#Compute eigenvalues on L as ndarray object
	L = nx.normalized_laplacian_matrix(graph)
	e = numpy.linalg.eigvals(L.A)
	#energy: squared sum of eigenvalues
	sum = 0
	for i in e:
		sum += i
	return_values.append(sum*sum)
	#Second largest eigenvalue
	
	#Number distict eigenvalues
	
	#Spectral Radius: largest |eigenvalue|
	return_values.append(max(e))
	return return_values
	
	
def labelAttributes(graph):
	return_values = []
	#Link Impurity: number of edges {u, v} where L(u) != L(v) / number of edges
	#Neighborhood Impurity: vertex degree for only those with different label
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
