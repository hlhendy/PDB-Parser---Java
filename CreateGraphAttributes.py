#inputs will be vector of atoms (nodes)
#and vector of contacts [0 or 1] (edges)
#will need to recreate loop to determine whether 
#the ith atom is connected to the jth atom

#or output tuples with atom1, atom2, contact

#read file, if 1 G.add_edge(atom1, atom2)
import networkx as nx
import csv
import sys
import time
import numpy.linalg
import Parser
import Distance
import GraphAttributes

##CONSTANTS
lRMSD_CRITERIA = 2
RES_DISTANCE = 4
BIN_CRITERIA = 8

#Attributes: 
#EIGENVALUE-BASED(0-3): energy, 2nd largest, #distict values, spectral radius
#LABELS-BASED(4-6): link impurity, neighborhood impurity, label entropy
#CLUSTER-BASED(7-10): closeness centrality, clustering coefficient, small-worldness criteria

def graphAttributes(graph):
	attributes = []
	attributes.append(GraphAttributes.eigenvalueAttributes(graph))
	attributes.append(GraphAttributes.labelAttributes(graph))
	attributes.append(GraphAttributes.clusterAttributes(graph))
	return attributes

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
	
	#output graph attributes for each data set
	dt = time.strftime("%Y%m%d-%H%M%S")
	with open('Output/'+output_prefix+dt+'.csv', 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerow(['energy', 'link_impurity', 'num_eigen', 'deg_imp', 'spectral_rad', 'label_entropy', 
			'closeness_centrality', 'clustering_coeff', 'second_eigen', 'small_worldness', 'near_native'])
		#Positive Data Set
		for i in range(len(withinlRMSD)):
			graph = nx.Graph()
			curr_conf = withinlRMSD[i]
			for j in range(len(curr_conf)-RES_DISTANCE):
				for k in range(j+RES_DISTANCE, len(curr_conf)):
					atom1 = curr_conf[j]
					atom2 = curr_conf[k]
					#add nodes to graph with labels
					graph.add_node(j, label=labels[j])
					graph.add_node(k, label=labels[k])
					#find euclidean distance between atoms
					d = Distance.euclideanDistance(atom1, atom2)
					#if less than BIN_CRITERIA, add edge
					if(d <= BIN_CRITERIA):
						graph.add_edge(j, k)
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
					graph.add_node(j, label=atom1[0])
					graph.add_node(k, label=atom2[0])
					#find euclidean distance between atoms
					d = Distance.euclideanDistance(atom1, atom2)
					#if less than BIN_CRITERIA, add edge
					if(d <= BIN_CRITERIA):
						graph.add_edge(j, k)
			#once graph is done, create attribute vector
			attributes = graphAttributes(graph)
			#add 0 to the end since decoy
			attributes.append(0)
			#and output to file as row
			writer.writerow(attributes)
		print("ATTRIBUTES HAVE BEEN OUTPUTTED")

if __name__ == "__main__":
   main(sys.argv[1:])
 
