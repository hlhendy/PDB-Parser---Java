##Calculate lRMSD and either output list of lRMSD values
##or return two lists sorted by withinLRMSD or moreThanLRMSD
##given a certain criteria
#import scipy.spatial
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
#import scipy as sp
from operator import itemgetter, attrgetter, methodcaller

##########################################################

#constants and global variables
RES_DISTANCE           = 5
start_time             = 0


##########################################################
### Distance Squared between two atoms
##########################################################
def distanceSquared(a1, a2):
	deltaX = a1[0] - a2[0]
	deltaY = a1[1] - a2[1]
	deltaZ = a1[2] - a2[2]
	result = (deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ)
	return result
  
def euclideanDistance(a1, a2):
    return math.sqrt(distanceSquared(a1, a2))

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
	atom: [aminoacid, x, y, z]
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

#Print distances from native
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

#Print pairwise distances
def printPairwiseDistances(distanceVectors, output_prefix):
	output_file = str(output_prefix) + '.distance'
	f = open(str(output_file), 'wt')
	for i in range(0, len(distanceVectors)):
		for j in range(i+1, len(distanceVectors)):
			print("I is" + i)
			distance = euclidean_distance(distanceVectors[i], distanceVectors[j])
			f.write(str(distance) + '\n')
	return

#Print lRMSDs from native 
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

#Sort conformations by lRMSDs value, according to criteria
def sortBylRMSDs(nativeconformation, conformations, criteria):
	withinlRMSD = []
	morethanlRMSD = []
	for i in range(0, len(conformations)):
		distance = lrmsd(nativeconformation[0], conformations[i])
		if distance <= criteria:
			withinlRMSD.append(conformations[i])
		else:
			morethanlRMSD.append(conformations[i])
	return withinlRMSD, morethanlRMSD




