"""

bead_analyze_csv_output_mask.py

Analyzes the raw quantification output from bead_mask_analyze_mod.py (or bead_collate_csv.py)
 - Sorts particles by condition
 - Performs background subtraction and error propagation across user-specified channels
 - Produces a human-readable spreadsheet output
 - If multiple channels, output normalized fluor ratios, assuming that the first listed 
   channel is the ratiometric control condition.


Modified for the internal masking and filtering protocol to not require sequentially numbered particles.

"""


import sys
import csv
from math import sqrt


def main(inCsv, channels):
	# Load the csv
	lines = []
	with open(inCsv, "r") as csvfile:
		reader = csv.reader(csvfile, delimiter=",")
		for row in reader:
			lines.append(row)
	#print(lines)

	# Identify all of the capture names.  They are start of index 0 up to first "_-" instance.
	captures = []
	for line in lines:
		if line[0][:5] == "Captu":
			index = line[0].find("_-")
			captures.append(line[0][:index])
	captures = list(set(captures))
	#print(captures)

	# Initialize a writer
	g = open(no_ext(inCsv)+"_analyzed.csv", "w")

	# Loop across the captures to perform analysis
	for capture in captures:
		#print(capture)
		# Find all the lines that correspond to this capture
		relevant_lines = []
		for i in range(0, len(lines)):
			if capture == lines[i][0][0:lines[i][0].find("_", lines[i][0].find("_")+1)]:
				relevant_lines.append(lines[i])

		# Instantiate a dictionary containing the data.
		capture_dict = {}
		for i in range(0, len(relevant_lines)):
			# Parse the micrograph name to get the number and channel info
			mcg_number, channel, prefix = up_to_the_channel(relevant_lines[i][0], channels)

			# Logic to figure out if this is a new micrograph or just a new channel
			if i == 0:
				new_mcg_flag = True
			elif relevant_lines[i] != relevant_lines[i-1]:
				new_mcg_flag = True
			else:
				new_mcg_flag = False

			# If this is a new micrograph, do the initial population of the dictionary
			if new_mcg_flag == True:
				if mcg_number not in capture_dict.keys():
					capture_dict[mcg_number] = {}
				
			# Is this the bg line?  If so, get the data and hold it.
			# It's going to get passed to the dictionary fill function along with the relevant lines.
			# This is essentially hacking the fact that the BG line is always written first in inCsv
			if "BG" in relevant_lines[i][1]:
				bg_mean, bg_stdev, bg_flag, bg_min, bg_max, bg_range = get_bg_vals(relevant_lines[i])
				capture_dict[mcg_number][channel] = {}
			else:
				particle_id = int(relevant_lines[i][1])
				capture_dict[mcg_number][channel][particle_id] = {}
				dict_pop(capture_dict[mcg_number][channel][particle_id], relevant_lines[i], mcg_number, channel, prefix, bg_mean, bg_stdev, bg_flag, bg_min, bg_max, bg_range)

		# Now that data is stored, it can be written out the way I want.
		# Capture-specific header and label info.
		g.write(capture+"\n")
		g.write("Parent Image, Particle ID, Mask Size, ")
		for j in range(0, len(channels)):
			g.write(channels[j]+" Mean, "+channels[j]+" Stdev, "+channels[j]+" Bg Mean, "+channels[j]+" Bg Stdev, "+channels[j]+" Corr. Mean, "+channels[j]+" Corr. Stdev, ")
			if j != 0:
				g.write(channels[j]+" norm to "+channels[0]+" mean, "+channels[j]+" norm err, ")
		g.write(channels[0]+" Bg error ratio, "+channels[0]+" Bg range\n")

		# Per-particle info
		for micrograph in capture_dict:
			#print(micrograph)
			# Define the range of particles to test
			try:
				highest_particle_id = max(list(capture_dict[micrograph][channels[0]].keys()))
				all_particle_ids = list(capture_dict[micrograph][channels[0]].keys())
			except:
				highest_particle_id = -1
			if highest_particle_id != -1:
				for i in range(0, len(all_particle_ids)):
					outStr = ""
					# Pass the info for every channel to the writer.
					outStr = outStr+str(micrograph)+", "+str(all_particle_ids[i])+", "+str(capture_dict[micrograph][channels[0]][all_particle_ids[i]]["mask_size"])+", "
					for j in range(0, len(channels)):
						#print(capture, micrograph, channels[j])
						outStr = outStr+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["mean"])+", "+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["stdev"])+", "+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["bg_mean"])+", "+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["bg_stdev"])+", "+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["mean_bg_corrected"])+", "+str(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["stdev_bg_corrected"])+", "
						# Calculate ratiometric values and add to the string
						if j != 0:
							mean_ratio = capture_dict[micrograph][channels[j]][all_particle_ids[i]]["mean_bg_corrected"] / capture_dict[micrograph][channels[0]][all_particle_ids[i]]["mean_bg_corrected"]
							stdev_ratio = propagate_err_div(capture_dict[micrograph][channels[j]][all_particle_ids[i]]["stdev_bg_corrected"], capture_dict[micrograph][channels[j]][all_particle_ids[i]]["mean_bg_corrected"], capture_dict[micrograph][channels[0]][all_particle_ids[i]]["stdev_bg_corrected"], capture_dict[micrograph][channels[0]][all_particle_ids[i]]["mean_bg_corrected"])
							outStr = outStr+str(mean_ratio)+", "+str(stdev_ratio)+", "
					outStr = outStr + str(capture_dict[micrograph][channels[0]][all_particle_ids[i]]["bg_flag"])+", "+ str(capture_dict[micrograph][channels[0]][all_particle_ids[i]]["bg_range"])+"\n"
					g.write(outStr)

	g.close()

	# Report back
	print("\nConversion complete, output file: "+no_ext(inCsv)+"_analyzed.csv")
 

def get_bg_vals(inLine):
	# Mean(3), Stdev(6), Flag(calc)
	mean = float(inLine[3])
	stdev = float(inLine[7])
	rel_var = stdev / mean
	min_val = int(float(inLine[5].strip()))
	max_val = int(float(inLine[6].strip()))
	range_val = max_val - min_val
	#print(min_val, max_val, range_val)
	return mean, stdev, rel_var, min_val, max_val, range_val


def dict_pop(inDict, inLine, mcgNumber, channel, prefix, bg_mean, bg_stdev, bg_flag, bg_min, bg_max, bg_range):
	# First, populate the gimmes - background values. Micrograph ID is an emergent property of the dict structure.
	inDict["bg_mean"] = float(bg_mean)
	inDict["bg_stdev"] = float(bg_stdev)
	inDict["bg_flag"] = bg_flag
	inDict["bg_min"] = bg_min
	inDict["bg_max"] = bg_max
	inDict["bg_range"] = bg_range

	# Now things we have to pull in manually
	inDict["mask_size"] = int(inLine[2])
	mean = float(inLine[3])
	inDict["mean"] = mean
	stdev = float(inLine[7])
	inDict["stdev"] = stdev

	# Finally, calculations
	inDict["mean_bg_corrected"] = mean - float(bg_mean)
	inDict["stdev_bg_corrected"] = propagate_err(stdev, float(bg_stdev))

	return inDict


def propagate_err_div(var1, mean1, var2, mean2):
	return sqrt((var1/mean1)**2 + (var2/mean2)**2)


def propagate_err(var1, var2):
	return sqrt(var1**2 + var2**2)


def up_to_the_channel(inStr, channels):
	search_key = ""
	for channel in channels:
		if channel in inStr:
			search_key = channel
			break
	if search_key == "":
		print("Could not parse this micrograph name: "+inStr)
		print("Exiting.")
		exit()
	else:
		#print(inStr)
		index = inStr.find(search_key)
		# Now need to find the number between the two underscores following the channel value.
		index2 = inStr.find("_", index)
		index3 = inStr.find(".", index2)
		mcg_number = int(inStr[index2+1:index3])
		prefix = inStr[:index]

		return mcg_number, search_key, prefix


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
	if len(sys.argv) >= 3:
		main(str(sys.argv[1]), sys.argv[2:])
	else:
		print("Check argument inputs: foo.py inCSV channel1 channel2 channel3 ...\n\n")
		exit()