"""
bead_analyze_everything_batches.py

This is a wrapper script to prep files for the bead_mask_analyze.py script.
 - Converts spaces in the filename into underscores.
 - Outputs a textfile called RUNME.sh that is the command to run to execute the python program.
 - The RUNME uses brute MPI technology to run batches of micrographs in separate python instances

"""


import sys
import os


def main(number_of_processes, channel, multichannel_flag):
	counter = 0
	names = []
	for root, dirs, files in os.walk("."):
		for filename in files:
			if filename[-5:] == ".tiff":
				if channel in filename:
					current_name = no_ext(filename)
					new_name = current_name.replace(" ", "_").replace("]", "_")
					os.rename(filename, new_name+".tiff")
					names.append(new_name+".tiff")
					counter += 1

	print("Processed "+str(counter)+" files.")

	# Define the number of micrographs to assign to each process
	#raw_per_process = float(counter) / float(number_of_processes)
	#max_per_process = int(raw_per_process) + 1

	# Split the filenames into batches of this size
	batches = []
	for i in range(0, number_of_processes):
		batches.append([])
	for i in range(1, len(names)+1):
		batches[(i % number_of_processes)].append(names[i-1])
	#print(batches)

	# Write out
	g = open("RUNME.sh", "w")
	for batch in batches:
		base_str = "python bead_mask_analyze_mod_mask.py "+multichannel_flag+" "
		for item in batch:
			base_str += item
			base_str += " "
		g.write(base_str + "& \n")
	g.write("\n")
	g.close()
				

def no_ext(inStr):
	"""
	Takes an input filename and returns a string with the file extension removed.

	"""
	prevPos = 0
	currentPos = 0
	while currentPos != -1:
		prevPos = currentPos
		currentPos = inStr.find(".", prevPos+1)
	return inStr[0:prevPos]


if __name__ == '__main__':	
	if len(sys.argv) == 4:
		main(int(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]))
	else:
		print("Check argument inputs: foo.py numProcesses channel multichannelFlag")
		exit()
