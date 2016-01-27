import networkx as nx
import csv
import sys
import time
import numpy.linalg
import Parser
import Distance
import GraphAttributes
from matplotlib import pyplot as plt

##CONSTANTS
lRMSD_CRITERIA = 3 #(between 2-4)
RES_DISTANCE = 4
BIN_CRITERIA = 8

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
	if len(argv) != 4:
		print('USAGE: <native pdb file> <pdb file> <model limit> <output file prefix>')
		sys.exit(2)
	try: #TODO: add better checking here
		native_in = str(argv[0])
		file_in = str(argv[1])
		nr_models = int(argv[2])
		output_prefix = str(argv[3])
	except:
		print('USAGE: <native pdb file> <pdb file> <model limit> <output file prefix>')
		sys.exit(2)
	#Create lists of conformations	
	labels, nativeconformation, conformations = Parser.PDB(native_in, file_in, nr_models)
	#Sort into positive and negative sets using lRMSD 
	withinlRMSD, morethanlRMSD = Distance.sortBylRMSDs(nativeconformation, conformations, lRMSD_CRITERIA)
	
	#output image of native graph
	nativeGraph = nx.Graph()
	curr_conf = nativeconformation[0]
	for j in range(len(curr_conf)-RES_DISTANCE):
		for k in range(j+RES_DISTANCE, len(curr_conf)):
			atom1 = curr_conf[j]
			atom2 = curr_conf[k]
			#add nodes to graph with labels
			nativeGraph.add_node(j)
			nativeGraph.node[j]['aminoAcid'] = labels[j]
			nativeGraph.add_node(k)
			nativeGraph.node[k]['aminoAcid'] = labels[k]
			#find euclidean distance between atoms
			d = Distance.euclideanDistance(atom1, atom2)
			#if less than BIN_CRITERIA, add edge
			if(d <= BIN_CRITERIA):
				nativeGraph.add_edge(j, k, distance=d)
	printGraph(nativeGraph, 'Output/PosGraphs/native')
	
	#output graph attributes for each data set
	#Note: removed newline='' from open() for linux
	dt = time.strftime("%Y%m%d-%H%M%S")
	with open('Output/'+output_prefix+dt+'.csv', 'w') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		writer.writerow(['num_edges', 'density','avg_degree','percent_endpoints','energy', 'second_eigen', 'unique_eigen', 'spectral_rad', 
			'link_impurity', 'neighborhood_impurity', 'avg_closeness', 'avg_clustering', 'small_worldness','eccentricity','diameter',
			'radius','%central_nodes', 'near_native'])
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
					graph.add_node(k)
					graph.node[k]['aminoAcid'] = labels[k]
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
					graph.add_node(k)
					graph.node[k]['aminoAcid'] = labels[k]
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
			#add 0 to the end since decoy
			attributes.append(0)
			#and output to file as row
			writer.writerow(attributes)
		print("ATTRIBUTES HAVE BEEN OUTPUTTED")

if __name__ == "__main__":
   main(sys.argv[1:])
 
