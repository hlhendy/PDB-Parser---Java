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
	native_atoms = Parser.countResidues(native_in)
	decoy_atoms = Parser.countResidues(file_in)
	if len(native_atoms) != len(decoy_atoms):
		print("Unequal, find longest common sequence")
		native_result, decoy_result = Parser.lcs(native_atoms, decoy_atoms)
		print("New length: " + str(len(native_result)))		
	else:
		native_result = []
		decoy_result = []	
	#Read and store native conformation
	nativelabels, nativeconformation = Parser.readConformations(str(native_in), 1, native_result)
	#Read decoys and store how many are within distance, morethan distance
	#using criteria{2,4}
	criteria = [2.0000000,4.0000000]
	models = 0
	atoms = []
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
				alpha_carbons = 1
			elif splt[0] == 'ATOM':
				if(splt[2] == 'CA'):
					if(len(decoy_result) == 0 or (str(splt[3]), int(splt[5])) in decoy_result):
						atoms.append((float(splt[6]), float(splt[7]), float(splt[8])))
			elif splt[0] == 'TER':
				if(len(atoms) > 0):
					currConf.append(atoms)
					models += 1
					distance = Distance.lrmsd(nativeconformation[0], currConf[0])
					#output_data.append([distance])
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
		#for d in output_data:	
		#	writer.writerow(d)
	print("Completed")
	

if __name__ == "__main__":
	main(sys.argv[1:])
