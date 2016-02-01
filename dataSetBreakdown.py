import Parser
import Distance

def main(argv):
	#Open files
	if len(argv) != 4:
		print('USAGE: <native pdb> <decoy pdb> <model limit> <output file prefix>.csv')
		sys.exit(2)
	try:
		native_in = str(argv[0])
		file_in = str(argv[1])
		nr_models = int(argv[2])
		output_prefix = str(argv[3])
	except:
		print('USAGE: <native pdb> <decoy pdb> <model limit> <output file prefix>.csv')
		sys.exit(2)
	#Count atoms and calculate LCS if needed
	native_atoms = countAtoms(native_in)
	decoy_atoms = countAtoms(file_in)
	if len(native_atoms) != len(decoy_atoms):
		native_result, decoy_result = lcs(native_atoms, decoy_atoms)
	else:
		native_result, decoy_result = []	
	#Read and store native conformation
	nativelabels, nativeconformation = readConformations(str(native_in), 1, native_result)
	#Read decoys and store how many are within distance, morethan distance
	#using criteria{2,4}
	criteria = [2,4]
	f_read = open(str(file_in), 'r')
	models = 0
	atoms = []
	nr_atoms = 0
	while models < nr_models:
		line = f_read.readline()
		splt = line.split()
		if splt[0] == 'MODEL':
			atoms = []
			nr_atoms = 0
		elif splt[0] == 'ATOM':
			nr_atoms += 1
			if(splt[2] == 'CA' and (len(lcs) == 0 or nr_atoms in lcs)):
				atoms.append((float(splt[6]), float(splt[7]), float(splt[])))
		elif splt[0] == 'TER':
			if(len(atoms) > 0):
				distance = Distance.lrmsd(nativeconformation, atoms)
				#add in code to increment some variables for each criteria#
	#Output results in table with protein name, lcs length, number within/morethan for each criteria
