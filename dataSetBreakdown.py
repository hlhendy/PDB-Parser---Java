import Parser
import Distance
import sys
import csv

def main(argv):
	if len(argv) != 3:
		print("Usage: <native.pdb> <decoys.pdb> <output.csv>")
		sys.exit(2)
	try:
		native_in = argv[0]
		file_in = argv[1]
		#nr_models = m
		output_file = argv[2]
	except:
		print("Usage: <native.pdb> <decoys.pdb> <output.csv>")
		sys.exit(2)
	#Count atoms and calculate LCS if needed
	native_atoms = Parser.countAtoms(native_in)
	decoy_atoms = Parser.countAtoms(file_in)
	if len(native_atoms) != len(decoy_atoms):
<<<<<<< HEAD
		print("Unequal, find longest common sequence")
=======
		print("Unequal: Native has " + str(len(native_atoms)) + " and conf has " + str(len(decoy_atoms)))
>>>>>>> 85092edccbc979a0ee1ca9e92795000d11b4cb86
		native_result, decoy_result = Parser.lcs(native_atoms, decoy_atoms)
	else:
		native_result = []
		decoy_result = []	
	#Read and store native conformation
	nativelabels, nativeconformation = Parser.readConformations(str(native_in), 1, native_result)
	#Read decoys and store how many are within distance, morethan distance
	#using criteria{2,4}
	criteria = [2,4]
	#f_read = open(str(file_in), 'r')
	models = 0
	atoms = []
	alpha_carbons = 1
	output_data = []
	currConf = []
	within2 = 0
	morethan2 = 0
	within4= 0
	morethan4 = 0
	with open(file_in, 'r') as f:
		for line in f:
		#while models < nr_models:
			#line = f_read.readline()
			splt = line.split()
			if splt[0] == 'MODEL':
				atoms = []
				currConf = []
<<<<<<< HEAD
				alpha_carbons = 1
			elif splt[0] == 'ATOM':
				if(splt[2] == 'CA'):
					if(len(decoy_result) == 0 or alpha_carbons in decoy_result):
						atoms.append((float(splt[6]), float(splt[7]), float(splt[8])))
					alpha_carbons += 1
=======
				nr_atoms = 0
			elif splt[0] == 'ATOM':
				nr_atoms += 1
				if(splt[2] == 'CA' and (len(decoy_result) == 0 or nr_atoms in decoy_result)):
					atoms.append((float(splt[6]), float(splt[7]), float(splt[8])))
>>>>>>> 85092edccbc979a0ee1ca9e92795000d11b4cb86
			elif splt[0] == 'TER':
				if(len(atoms) > 0):
					currConf.append(atoms)
					models += 1
					distance = Distance.lrmsd(nativeconformation[0], currConf[0])
					if distance <= criteria[0]:
						within2 += 1
						within4 += 1
					else:
						morethan2 += 1
						if distance <= criteria[1]:
							within4 +=1 
						else: 
							morethan4 += 1
		#Output results in table with protein name, lcs length, number within/morethan for each criteria
		output_data.append(native_in[5:-4])
		output_data.append(len(nativeconformation[0]))
		output_data.append(within2+morethan2)
		output_data.append(within2)
		output_data.append(morethan2)
		output_data.append(within4)
		output_data.append(morethan4)
	with open(output_file, 'a+') as csvfile:
		writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
		if(csvfile.readline() == ""):
			writer.writerow(["Protein", "Num CA", "Num Confs", "Within 2", "Morethan 2", "Within 4", "Morethan 4"])	
		writer.writerow(output_data)
	print("Completed")
	
#def main(argv):
#	if len(argv) != 1:
#		print("Please enter name of dat file")
#		sys.exit(2)
#	try:
#		dat_file = argv[0]
#	except:
#		print("Please enter name of dat file")
#		sys.exit(2)
#	#Data file should include number of proteins as first line
#	#Each line should contain <native file> <conf file> <number of models>
#	output_file = 'Output/DatSetBreakdown.csv'
#	with open(output_file, 'a') as csvfile:
#		writer = csv.writer(csvfile, delimiter=',', quoting=csv.QUOTE_MINIMAL)
#		writer.writerow(['Protein Name', 'LCS', 'Within 2', 'Morethan 2', 'Within 4', 'Morethan 4'])
#	#For each line, assign argument variables and send to dataSetBreakdown
#	#while(next line):
#		#line = line.next
#		#native_in = split(line)[0] ... and so on
#		dataSetBreakdown.breakdown(native_in, file_in, nr_models, output_file)

if __name__ == "__main__":
	main(sys.argv[1:])
