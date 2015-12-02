##Modified for Python3  - Heather Hendy OCTOBER 2015
import scipy.spatial
import math
import sys, getopt
import random
import collections
from collections import defaultdict
import copy
import numpy as np
from numpy.linalg import det, svd
import os
import time
import scipy as sp
from operator import itemgetter, attrgetter, methodcaller

##########################################################

#constants and global variables
RES_DISTANCE           = 5
start_time             = 0

###################################################################
## read file with models and file with corresponding energies
## reads all atoms and coordinates
## pass conformations/models and energies to graph building
###################################################################
def readConformations(file_in, nr_models):
	f_read = open(str(file_in), 'r')
	models = 0
	conformations = []
	atoms = []
	atom_chars = []
    #read models
	while models < nr_models:
		line = f_read.readline()
		splt = line.split()
		#parse lines
		if splt[0] == 'MODEL':
			atoms = []
		elif splt[0] == 'ATOM': #and (splt[2] == 'CA'):
			atoms.append((str(splt[2]),float(splt[6]), float(splt[7]), float(splt[8])))
		elif splt[0] == 'TER': # or splt[0] == 'END' or splt[0] == 'ENDMDL':
			if(len(atoms) > 0):
				conformations.append(atoms)
				models +=1
			else:
				print("0 atoms at model ", models)
			if models % (models * .10) == 0:
				print("read model", models, "with ", len(atoms), "CA atoms")
	return conformations

#############################################################################
### Calculate Euclidean distance between two vectors
# this is minkowski norm 2, which we use in the kd tree here
#############################################################################
def euclidean_distance(x, y):
	if(len(x) != len(y)):
		print("at euclidean: ", len(x), "!=", len(y), "!!")
		sys.exit()
	d = 0
	for i in range(0, len(x)):
		d = d + (x[i] - y[i])**2
	return math.sqrt(d)

#this is minkowski norm 1
def manhattan_distance(x, y):
	if(len(x) != len(y)):
		sys.exit()
	d = 0
	for i in range(0, len(x)):
		d = d + math.abs(x[i] - y[i])
	return d

###################################################################
## build distance vectors
###################################################################
def buildDistanceVectorsFromConfs(conformations):
	nr_confs = len(conformations)
	nr_atoms = len(conformations[0])

	d_vectors = []
	for model in range(0, nr_confs):
		if(len(conformations[model]) < nr_atoms):
			print(model, len(conformations[model]))
		distance_vector = []
		k = 0
		for i in range(0, nr_atoms-RES_DISTANCE):
			for j in range(i+RES_DISTANCE, nr_atoms):
				distance_vector.append([])
				distance_vector[k] = euclidean_distance(conformations[model][i], conformations[model][j])
				k += 1
		dv = np.array(distance_vector)
		d_vectors.append(np.array(distance_vector))

		if(model % 1000 == 0):
			print("built distance vector for conformation ", model)
		dv2 = np.array(d_vectors)
	return dv2


def centroids(conformation1, conformation2, nr_atoms):
	"""
	centers a conformation at the centroid
	"""
	center1 = [0, 0, 0]
	center2 = [0, 0, 0]
	for atom in range(0, nr_atoms):
		center1[0] += conformation1[atom][0]
		center1[1] += conformation1[atom][1]
		center1[2] += conformation1[atom][2]

		center2[0] += conformation2[atom][0]
		center2[1] += conformation2[atom][1]
		center2[2] += conformation2[atom][2]
	for ww in range(0,3):
		center1[ww] /= nr_atoms
		center2[ww] /= nr_atoms
	return center1, center2


############################################################################
### conformations moved so their centroids are at (0,0,0)
############################################################################
def realign(center1, center2, conformation1, conformation2, nr_atoms):
	n_conf1 = []
	n_conf2 = []
	app1 = n_conf1.append
	app2 = n_conf2.append
	for atom in range(0, nr_atoms):
		app1(list(conformation1[atom]))
		app2(list(conformation2[atom]))
		for x in range(0, 3):
			n_conf1[atom][x] -= center1[x]
			n_conf2[atom][x] -= center2[x]
	return n_conf1, n_conf2
	

#############################################################################
### one should be chosen as reference
### assumption: first argument used as reference
#############################################################################
def lrmsd(conf1, conf2):
	"""
	Calculates lrmsd between two instances of two proteins (keys from create_conformations output dictionary)
	
	Arguments:
	conf1 -- first protein conformation (tuple of tuples) -- turn into list of lists
	conf2 -- second protein conformation (tuples of tuples) -- turn into list of lists
	
	Returns:
	number representing the lrmsd of the two conformations
	"""

	nr_atoms = len(conf1)
	center1, center2 = centroids(conf1, conf2, nr_atoms)

	#Move conf1 and conf2 so their centroids are at origin
        #results stored in conf1_prime and conf2_prime

	conf1_prime, conf2_prime = realign(center1, center2, conf1, conf2, nr_atoms)

	#s0 = time.time()	

        #calculate optimal rotation
	X = np.array(conf1_prime).T
	Y = np.array(conf2_prime).T
	M = np.dot(X, Y.T)
	V, S, Wt = svd(M)
	if det(M) > 0:
		U = np.dot(Wt.T, V.T)
	else:
		new = np.array([[1,0,0],[0,1,0],[0,0,-1]])
		U_p = np.dot(Wt.T, new)
		U = np.dot(U_p, V.T)
	xp = np.dot(U, X)
	sc = np.sum((xp-Y)**2)
	sc = math.sqrt(sc/nr_atoms)
	return sc

def printDistancesFromNative(nativedistanceVector, distanceVectors, output_prefix):
	distancesFromNative = []
	output_file = str(output_prefix) + '.distance'
	f = open(str(output_file), 'wb')
	print("opened ", str(output_file))
	for i in range(0, len(distanceVectors)):
		print(i, len(distanceVectors))
		distance = euclidean_distance(nativedistanceVector[0], distanceVectors[i])
		print(i, distance)
		f.write(str(distance) + '\n')
		distancesFromNative.append(distance)
	return distancesFromNative


def printPairwiseDistances(distanceVectors, output_prefix):
	output_file = str(output_prefix) + '.distance'
	f = open(str(output_file), 'wb')
	for i in range(0, len(distanceVectors)):
		for j in range(i+1, len(distanceVectors)):
			print("I is" + i)
			distance = euclidean_distance(distanceVectors[i], distanceVectors[j])
			f.write(str(distance) + '\n')
	return

def printLRMSDsFromNative(nativeconformation, conformations, output_prefix):
	lrmsdsFromNative = []
	output_file = str(output_prefix) + '.rmsd'
	f = open(str(output_file), 'wb')
	for i in range(0, len(conformations)):
		distance = lrmsd(nativeconformation[0], conformations[i])
		f.write(str(distance) + '\n')
	if i % 100 == 0:
		print("processed lrmsds of conformation ", i)
	lrmsdsFromNative.append(distance)
	return lrmsdsFromNative

def printPairwiseLRMSDs(conformations, output_prefix):
	output_file = str(output_prefix) + '.rmsd'
	f = open(str(output_file), 'wb')
	for i in range(0, len(conformations)):
		for j in range(i+1, len(conformations)):
			distance = lrmsd(conformations[i], conformations[j])
			f.write(str(distance) + '\n')
		if i % 100 == 0:
			print("processed lrmsds of conformation ", i)
	return

#accept native and decoy conformations
#find longest common subsequence
#return arrays of indices indicating locations of lcs elements in each array, a and b
def lcs(a, b):
	C = [[0 for j in range(len(b)+1)] for i in range(len(a)+1)] 
	# row 0 and column 0 are initialized to 0 already
	for i, x in enumerate(a):
		for j, y in enumerate(b):
			if x == y:
				C[i+1][j+1] = C[i][j] + 1
			else:
				C[i+1][j+1] = max(C[i+1][j], C[i][j+1])
	# read the substring out from the matrix
	result_a = [len(C)]
	result_b = [len(C)]
	x, y = len(a), len(b)
	r = 0
	while x != 0 and y != 0:
		if C[x][y] == C[x-1][y]:
			x -= 1
		elif C[x][y] == C[x][y-1]:
			y -= 1
		else:
			assert a[x-1] == b[y-1]
			result_a[r] = x
			result_b[r] = y
			r += 1
			x -= 1
			y -= 1
	return result_a, result_b
	
def extractCAs(list):
	#iterate through list of conformations
	print(len(list[0]))
	for i in range(0,len(list)):
		#iterate through list of atoms in each conformation
		for j in range(0, len(list[i])):
		#remove any atom that does not have CA as name
			if list[i][j][0] != 'CA':
				del list[i]
	#return list
	print(len(list[0]))
	return list

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
	global start_time
	start_time = time.time()
	print("READING INPUT FILES")
	nativeconformation = readConformations(str(native_in), 1)
	print("READ NATIVE CONF OF ", len(nativeconformation[0]), "CA atoms")
	conformations = readConformations(str(file_in), nr_models)
	print("READ ", len(conformations), " CONFORMATIONS (", len(conformations[0]), "CA atoms each) -- time passed:", time.time() - start_time, "(s)\n")
	##############################################################################

	##############################################################################
	#determine longest common subsequence if needed
	print("DETERMINING LONGEST COMMON SUBSEQUENCE")
	if len(nativeconformation[0]) != len(conformations[0]):
		print("Lengths are different: " + str(len(nativeconformation[0])) + " " + str(len(conformations[0])))
		nativeconformation, conformations = lcs(nativeconformation, conformations)
	##############################################################################
	
	##############################################################################
	#finding alpha carbons
	print("REDUCING TO ALPHA CARBONS")
	nativeconformation = extractCAs(nativeconformation)
	conformations = extractCAs(conformations)
	##############################################################################

	##############################################################################
	#reading building distance vectors
	print("BUILDING DISTANCE REPRESENTATIONS OF CONFORMATIONS")
	nativedistanceVector = buildDistanceVectorsFromConfs(nativeconformation)
	print("BUILT ", len(nativedistanceVector), " DISTANCE VECTORS (", len(nativedistanceVector[0]), "distances each) -- time passed:", time.time() - start_time, "(s)\n")
	distanceVectors = buildDistanceVectorsFromConfs(conformations)
	print("BUILT ", len(distanceVectors), " DISTANCE VECTORS (", len(distanceVectors[0]), "distances each) -- time passed:", time.time() - start_time, "(s)\n")
	##############################################################################

	##############################################################################
	#print distances
	print("PRINTING DISTANCES FROM NATIVE")
	distancesFromNative = printDistancesFromNative(nativedistanceVector, distanceVectors, output_prefix)
	print("PRINTED DISTANCES -- time passed:", time.time() - start_time, "(s)\n")
	##############################################################################

	##############################################################################
	#print lrmsds
	print("PRINTING LRMSDS FROM NATIVE")
	lrmsdsFromNative = printLRMSDsFromNative(nativeconformation, conformations, output_prefix)
	print("PRINTED LRMSDs -- time passed:", time.time() - start_time, "(s)\n")
	##############################################################################

if __name__ == "__main__":
   main(sys.argv[1:])


#b_c_e_d("1PGB_RosettaDecoys.pdb", 1098, "1PGB_RosettaDecoys.score12", "1PGB_Results.txt", 3.0)


