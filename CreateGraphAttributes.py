import networkx as nx
import csv
import sys
import time
import numpy.linalg
import Parser
import Distance
import GraphAttributes
#from matplotlib import pyplot as plt //In python 2.6.6 but not in python 3.3

##CONSTANTS
#lRMSD_CRITERIA = 4 #produces more balanced pos/neg data points than 2
RES_DISTANCE = 2 #Largest distance where all graphs are connected
BIN_CRITERIA = 8 #In order to put an edge between nodes (u, v), the distance between them must be 8A or less
HPHOBIC = ['ALA', 'ILE', 'LEU', 'PHE', 'VAL', 'PRO', 'GLY']
HPHILIC = ['ARG', 'LYS', 'ASP', 'GLU', 'GLN', 'ASN', 'HIS', 
			'SER', 'THR', 'TYR', 'CYS', 'MET', 'TRP']

#Attributes(18):
#GENERAL: 		number of edges, avg. degree, density, %endpoints 
#EIGENVALUE-BASED: 	energy, 2nd largest, #distict values, spectral radius
#LABELS-BASED: 		link impurity, neighborhood impurity, label entropy
#CLUSTER-BASED: 	closeness centrality, clustering coefficient
#SMALL-WORLDNESS:	small-worldness criteria
#ECCENTRICITY-BASED:	avg eccentricity, diameter, radius, %central nodes
def graphAttributes(graph):
	attributes = []
	for a in GraphAttributes.generalAttributes(graph):
		attributes.append(a)
	for a in GraphAttributes.eigenvalueAttributes(graph):
		attributes.append(a)
	for a in GraphAttributes.labelAttributes(graph):
		attributes.append(a)
	for a in GraphAttributes.clusterAttributes(graph):
		attributes.append(a)
	for a in GraphAttributes.smallWorldness(graph):
		attributes.append(a)
	for a in GraphAttributes.eccentricityAttributes(graph):
		attributes.append(a)
	return attributes
	
def printGraph(graph, filename):
	G = graph
	pos = nx.spring_layout(G)
	nx.draw(G, pos)
	node_labels = nx.get_node_attributes(G,'aminoAcid')
	nx.draw_networkx_labels(G, pos, labels = node_labels)
	#edge_labels = nx.get_edge_attributes(G, 'distance')
	#nx.draw_networkx_edge_labels(G, pos, labels = edge_labels)
	plt.savefig(filename + '.png')
	#plt.show()

#######################################################################################
def main(argv):
	if len(argv) != 5:
		print('USAGE: <native pdb file> <pdb file> <model limit> <output file prefix> <lrmsd criteria>')
		sys.exit(2)
	try: #TODO: add better checking here
		native_in = str(argv[0])
		file_in = str(argv[1])
		nr_models = int(argv[2])
		output_prefix = str(argv[3])
		lrmsd_criteria = int(argv[4])
	except:
		print('USAGE: <native pdb file> <pdb file> <model limit> <output file prefix> <lrmsd criteria>')
		sys.exit(2)
	#Create lists of conformations	
	labels, nativeconformation, conformations = Parser.PDB(native_in, file_in, nr_models)
	#Sort into positive and negative sets using lRMSD 
	withinlRMSD, morethanlRMSD = Distance.sortBylRMSDs(nativeconformation, conformations, lrmsd_criteria)
	
	#output image of native graph
	#nativeGraph = nx.Graph()
	#curr_conf = nativeconformation[0]
	#for j in range(len(curr_conf)-RES_DISTANCE):
	#	for k in range(j+RES_DISTANCE, len(curr_conf)):
	#		atom1 = curr_conf[j]
	#		atom2 = curr_conf[k]
	#		#add nodes to graph with labels
	#		nativeGraph.add_node(j)
	#		nativeGraph.node[j]['aminoAcid'] = labels[j]
	#		nativeGraph.add_node(k)
	#		nativeGraph.node[k]['aminoAcid'] = labels[k]
	#		#find euclidean distance between atoms
	#		d = Distance.euclideanDistance(atom1, atom2)
	#		#if less than BIN_CRITERIA, add edge
	#		if(d <= BIN_CRITERIA):
	#			nativeGraph.add_edge(j, k, distance=d)
	#printGraph(nativeGraph, 'Output/PosGraphs/native')
	
	#output graph attributes for each data set
	#Note: removed newline='' from open() for linux
	dt = time.strftime("%Y%m%d-%H%M%S")
	with open('Output/'+output_prefix+dt+'.csv', 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		writer.writerow(['num_edges', 'density','avg_degree','percent_endpoints','energy', 'second_eigen', 'unique_eigen', 'spectral_rad', 'inverse_product',
			'link_impurity', 'neighborhood_impurity', 'avg_closeness', 'avg_clustering', 'small_worldness','eccentricity','diameter',
			'radius','%central_nodes', '%Hydrophobic_center', 'near_native'])
		#Positive Data Set
		for i in range(len(withinlRMSD)):
			graph = nx.Graph()
			curr_conf = withinlRMSD[i]
			for j in range(len(curr_conf)-RES_DISTANCE):
				for k in range(j+RES_DISTANCE, len(curr_conf)):
					atom1 = curr_conf[j]
					atom2 = curr_conf[k]
					#add nodes to graph with labels
					graph.add_node(j)
					graph.node[j]['aminoAcid'] = labels[j]
					if(labels[j] in HPHOBIC):
						graph.node[j]['hydro'] = 'phobic'
					else:
						graph.node[j]['hydro'] = 'philic'
					graph.add_node(k)
					graph.node[k]['aminoAcid'] = labels[k]
					if(labels[j] in HPHOBIC):
						graph.node[k]['hydro'] = 'phobic'
					else:
						graph.node[k]['hydro'] = 'philic'
					#find euclidean distance between atoms
					d = Distance.euclideanDistance(atom1, atom2)
					#if less than BIN_CRITERIA, add edge
					if(d <= BIN_CRITERIA):
						graph.add_edge(j, k, distance=d)
			##FOR TESTING ONLY
			#printGraph(graph, 'Output/PosGraphs/pos_'+str(i))
			#################
			#once graph is done, create attribute vector
			attributes = graphAttributes(graph)
			##FOR TESTING##
			#attributes = []
			#if(not nx.is_connected(graph)):
			#	print("Graph " + i + "from within is not connected")
			#	sys.exit(2)
			#else:
			#	attributes.append(nx.is_connected(graph))
			#add 1 to the end since near native
			attributes.append(1)
			#and output to file as row
			writer.writerow(attributes)
		#Negative Data Set
		for i in range(len(morethanlRMSD)):
			graph = nx.Graph()
			curr_conf = morethanlRMSD[i]
			for j in range(len(curr_conf)-RES_DISTANCE):
				for k in range(j+RES_DISTANCE, len(curr_conf)):
					atom1 = curr_conf[j]
					atom2 = curr_conf[k]
					#add nodes to graph with labels
					graph.add_node(j)
					graph.node[j]['aminoAcid'] = labels[j]
					if(labels[j] in HPHOBIC):
						graph.node[j]['hydro'] = 'phobic'
					else:
						graph.node[j]['hydro'] = 'philic'
					graph.add_node(k)
					graph.node[k]['aminoAcid'] = labels[k]
					if(labels[k] in HPHOBIC):
						graph.node[k]['hydro'] = 'phobic'
					else:
						graph.node[k]['hydro'] = 'philic'
					#find euclidean distance between atoms
					d = Distance.euclideanDistance(atom1, atom2)
					#if less than BIN_CRITERIA, add edge
					if(d <= BIN_CRITERIA):
						graph.add_edge(j, k, distance=d)
			##FOR TESTING ONLY
			#printGraph(graph, 'Output/NegGraphs/neg_'+str(i))
			#################
			#once graph is done, create attribute vector
			attributes = graphAttributes(graph)
			##FOR TESTING ONLY##
			#if(not nx.is_connected(graph)):
			#	print("Graph " + i + "from morethan is not connected")
			#	sys.exit(2)
			#else:
			#	attributes.append(nx.is_connected(graph))
			#add 0 to the end since decoy
			attributes.append(0)
			#and output to file as row
			writer.writerow(attributes)
		print("ATTRIBUTES HAVE BEEN OUTPUTTED")

if __name__ == "__main__":
   main(sys.argv[1:])
 
