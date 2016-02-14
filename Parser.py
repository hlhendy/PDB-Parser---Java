######################################################################
##	Parser												
##		- PDB(String filename): reads PDB file and outputs lists of
##			atoms (tuples that include atom name, amino acid name, 
##			coords)
######################################################################
import time
import numpy as np
from collections import Counter
#constants and global variables
START_TIME = 0

######################################################################
## read in atom names and return list of chars for each atom name
######################################################################
def countResidues(file_in):
        f_read = open(str(file_in), 'r')
        res_names = []
        models = 0
        #read first model and store atom names only
        while models < 1:
                line = f_read.readline()
                splt = line.split()
                #parse lines
                if splt[0] == 'ATOM' and splt[2] == 'CA': #changed to capture only CA atoms
                        res_names.append((str(splt[3]), int(splt[5]))) #Changed to capture residue name/number rather than atom name
                elif splt[0] == 'TER': # or splt[0] == 'END' or splt[0] == 'ENDMDL':
                        if len(res_names) > 0:
                                models +=1
                        else:
                                print("0 atoms at model ", models)
        f_read.close()
	print("Length is: " + str(len(res_names)))
        return res_names

###################################################################
## read file with models and file with corresponding energies
## reads all atoms and coordinates, stores alpha carbons only (if in LCS)
## pass conformations/models and energies to graph building
###################################################################
def readConformations(file_in, nr_models, lcs):
	f_read = open(str(file_in), 'r')
	models = 0
	alpha_carbons = 1
	conformations = []
	atoms = []
	labels = []
	print("Starting New Conf")
   #Parse and store only CA atoms
	while models < nr_models:
		line = f_read.readline()
		splt = line.split()
		#parse lines
		if(splt[0] == 'MODEL'):
			atoms = []
			alpha_carbons = 1
		elif splt[0] == 'ATOM':
			if(splt[2] == 'CA'):
				if(len(lcs) == 0 or (str(splt[3]), int(splt[5])) in lcs):
					#Store coords for atom
					atoms.append((float(splt[6]), float(splt[7]), float(splt[8])))
					#store amino acid names as labels
					#only need to do this once
					if(models < 1):
						labels.append(str(splt[3]))
				alpha_carbons += 1
		elif splt[0] == 'TER': # or splt[0] == 'END' or splt[0] == 'ENDMDL':
			if(len(atoms) > 0):
				conformations.append(atoms)
				models +=1
			else:
				print("0 atoms at model ", models)
			if models % 10000 == 0:
				print("read model", models, "with ", len(atoms), "CA atoms")
	return labels, conformations

#accepts list of CAs from native and first decoy (tuples with res name,number)
#find longest common subsequence
#return LCS with tuples names/number
def lcs(X , Y):
    # find the length of the strings
	m = len(X)
	n = len(Y)
    # declaring the arrays with zeros
	C = np.zeros((m+1, n+1))
	dir = np.zeros((m+1, n+1))
	for i in range(m+1):
		for j in range(n+1):
			if X[i-1][0] == Y[j-1][0]: #Removed Counter
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
	while(i>0 and j>0):
		if X[i-1][0] == Y[j-1][0]: #Removed Counter
			LCS = str(X[i-1]) + "-" + LCS
			#Return index of residue name
			result_X.append(X[i-1])
			result_Y.append(Y[j-1])
			i-= 1
			j-= 1
			length -= 1
		elif dir[i][j] == 0:
			#LCS[C[i][j]] = X[i]
			i-= 1
		else:
			#lCS[C[i][j]] = Y[j]
			j-= 1	
	return result_X, result_Y

#######################################################################################
def PDB(native_in, file_in, nr_models):
	global START_TIME
	START_TIME = time.time()
	print("DETERMINING LENGTH OF ATOMS")
	native_atoms = countResidues(native_in)
	decoy_atoms = countResidues(file_in)
	#determine longest common subsequence if needed
	if len(native_atoms) != len(decoy_atoms):
		print("Lengths are different: " + str(len(native_atoms)) + " " + str(len(decoy_atoms)))
		print("DETERMINING LONGEST COMMON SUBSEQUENCE")
		native_result, decoy_result = lcs(native_atoms, decoy_atoms)
		print("Lengths are now: " + str(len(native_result)) + " " + str(len(decoy_result)))
	else:
                print("Lengths are the same: " + str(len(native_atoms)) + " " + str(len(decoy_atoms)))
		native_result = []
		decoy_result = []

	#send readConformations array of indices from LCS if applicable, send empty list if not needed
	#use this list to include only CAs that are in the LCS
	print("READING INPUT FILES")
	nativelabels, nativeconformation = readConformations(str(native_in), 1, native_result)
	print("READ NATIVE CONF OF ", len(nativeconformation[0]), "CA atoms")
	labels, conformations = readConformations(str(file_in), nr_models, decoy_result)
	if(nativelabels != labels):
		#do something
		print("labels don't match!!")
	else:
		labels_final = labels
	print("READ ", len(conformations), " CONFORMATIONS (", len(conformations[0]), "CA atoms each) -- time passed:", time.time() - START_TIME, "(s)\n")
	return labels, nativeconformation, conformations;
	##############################################################################

	



