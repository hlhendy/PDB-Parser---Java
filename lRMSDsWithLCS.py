##Modified for Python3  - Heather Hendy OCTOBER 2015
##Outputs lRMSDs for use in labeling positive and negative set
##Added ability to calculate Longest Common Subsequence
import scipy.spatial
import math
import sys, getopt
import random
import collections
from collections import defaultdict
from collections import Counter
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
## read in atom names and return list of chars for each atom name
###################################################################
def countAtoms(file_in):
        f_read = open(str(file_in), 'r')
        atom_chars = []
        models = 0
        #read first model and store atom names only
        #TODO: Strip to return only char for name
        while models < 1:
                line = f_read.readline()
                splt = line.split()
                #parse lines
                if splt[0] == 'ATOM':
                        atom_chars.append(str(splt[2]))
                elif splt[0] == 'TER': # or splt[0] == 'END' or splt[0] == 'ENDMDL':
                        if len(atom_chars) > 0:
                                models +=1
                        else:
                                print("0 atoms at model ", models)
        f_read.close()
        return atom_chars

###################################################################
## read file with models and file with corresponding energies
## reads all atoms and coordinates, stores alpha carbons only (if in LCS)
## pass conformations/models and energies to graph building
###################################################################
def readConformations(file_in, nr_models, lcs):
	f_read = open(str(file_in), 'r')
	models = 0
	nr_atoms = 0
	conformations = []
	atoms = []
	print("Starting New Conf")
   #Parse and store only CA atoms
	while models < nr_models:
                line = f_read.readline()
                splt = line.split()
                #parse lines
                if splt[0] == 'MODEL':
                        atoms = []
                        nr_atoms = 0
                elif splt[0] == 'ATOM':
                        nr_atoms += 1
                        if splt[2] == 'CA' and (len(lcs) == 0 or nr_atoms in lcs):
                                #print(str(nr_atoms) + " " + str(splt[2]) + "\n")
                                atoms.append((float(splt[6]), float(splt[7]), float(splt[8])))
                elif splt[0] == 'TER': # or splt[0] == 'END' or splt[0] == 'ENDMDL':
                        if(len(atoms) > 0):
                                conformations.append(atoms)
                                models +=1
                        else:
                                print("0 atoms at model ", models)
                        if models % 100 == 0:
                                print("read model", models, "with ", len(atoms), "CA atoms")

	#f_read.close()
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
	f = open(output_file, 'wt')
	print("opened ", str(output_file))
	for i in range(0, len(distanceVectors)):
		print(i, len(distanceVectors))
		distance = euclidean_distance(nativedistanceVector[0], distanceVectors[i])
		print(i, distance)
		f.write(str(distance)+'\n')
		distancesFromNative.append(distance)
	return distancesFromNative


def printPairwiseDistances(distanceVectors, output_prefix):
	output_file = str(output_prefix) + '.distance'
	f = open(str(output_file), 'wt')
	for i in range(0, len(distanceVectors)):
		for j in range(i+1, len(distanceVectors)):
			print("I is" + i)
			distance = euclidean_distance(distanceVectors[i], distanceVectors[j])
			f.write(str(distance) + '\n')
	return

def printLRMSDsFromNative(nativeconformation, conformations, output_prefix):
	lrmsdsFromNative = []
	output_file = str(output_prefix) + '.rmsd'
	f = open(output_file, 'wt')
	for i in range(0, len(conformations)):
		distance = lrmsd(nativeconformation[0], conformations[i])
		f.write(str(distance) + '\n')
	if i % 100 == 0:
		print("processed lrmsds of conformation ", i)
	lrmsdsFromNative.append(distance)
	return lrmsdsFromNative

def printPairwiseLRMSDs(conformations, output_prefix):
	output_file = str(output_prefix) + '.rmsd'
	f = open(str(output_file), 'wt')
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
def lcs(X , Y):
    # find the length of the strings
	m = len(X)
	n = len(Y)
    # declaring the arrays with zeros
	C = np.zeros((m+1, n+1))
	dir = np.zeros((m+1, n+1))
	for i in range(m+1):
		for j in range(n+1):
			if Counter(X[i-1]) == Counter(Y[j-1]):
				C[i][j] = C[i-1][j-1]+1
				dir[i][j] = 1 #diagonal
			else:
				C[i][j] = max(C[i-1,j], C[i,j-1])
				if C[i][j] == C[i-1][j]:
					dir[i][j] = 0 #up
				else:
					dir[i][j] = -1 #down

    # C[m][n] contains the length of LCS
	length = C[m][n]
	result_X = [] #indices for X
	result_Y = [] #indices for Y
	LCS = ""
	i = m
	j = n
	while(i>0 and j >0):
		if Counter(X[i-1]) == Counter(Y[j-1]):
			LCS = X[i-1] + "-" + LCS
			result_X.append(i)
			result_Y.append(j)
			i-= 1
			j-= 1
			length -= 1
		elif dir[i][j] == 0:
			#LCS[C[i][j]] = X[i]
			i-= 1
		else:
			#lCS[C[i][j]] = Y[j]
			j-= 1
##	file = open('LCS_output', 'w')
##	file.write(str(LCS) + "\n")
##	file.write(str(result_X))
	return result_X, result_Y

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
	print("DETERMINING LENGTH OF ATOMS")
	native_atoms = countAtoms(native_in)
	decoy_atoms = countAtoms(file_in)
	#determine longest common subsequence if needed
	if len(native_atoms) != len(decoy_atoms):
		print("Lengths are different: " + str(len(native_atoms)) + " " + str(len(decoy_atoms)))
		print("DETERMINING LONGEST COMMON SUBSEQUENCE")
		native_result, decoy_result = lcs(native_atoms, decoy_atoms)
		print("Lengths are now: " + str(len(native_result)) + " " + str(len(decoy_result)))
	else:
                native_result, decoy_result = []

	#send readConformations array of indices from LCS if applicable, send empty list if not needed
	#use this list to include only CAs that are in the LCS
	print("READING INPUT FILES")
	nativeconformation = readConformations(str(native_in), 1, native_result)
	print("READ NATIVE CONF OF ", len(nativeconformation[0]), "CA atoms")
	conformations = readConformations(str(file_in), nr_models, decoy_result)
	print("READ ", len(conformations), " CONFORMATIONS (", len(conformations[0]), "CA atoms each) -- time passed:", time.time() - start_time, "(s)\n")
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


