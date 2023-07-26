"""

bead_calc_summary_stats.py

Take the filtered csv file and calculate summary statistics for each condition.

Name		Mean Norm. Signal	Mean Norm. Singal. Err (STDEV.S)	N
Capture_7
Capture_8
Capture_9

Replicate groups is anything in csv format that lists comma-separated groups to be merged on each line

"""


import sys
import numpy as np


def main(in_csv, replicate_groups):
	# Load the csv
	with open(in_csv, "r") as csvfile:
		lines = csvfile.readlines()

	# Identify all the captures and their indeces
	capture_indeces = []
	capture_names = []
	for i in range(0, len(lines)):
		items = lines[i].split(",")
		if items[0][:5] == "Captu":
			capture_indeces.append(i)
			capture_names.append(items[0].replace("\n", ""))

	# For each capture, walk forward per particle until you find another line len < 4
	
	vals_dict = {}
	for j in range(0, len(capture_indeces)):
		# From each particle, grab the normalized signal value at index 15 and put it in a list
		norm_signals = []
		for i in range(capture_indeces[j]+2, len(lines)):
			items = lines[i].split(",")
			if len(items) < 4 or items[0] == "\n" or "Parent Image" in items[0] or "Capture" in items[0]:
				break
			#print(lines[i])
			#print(i, lines[i][0], lines[i][1], lines[i][2])
			norm_signals.append(float(items[15]))
		
		if len(norm_signals) != 0:
			# Load a dict with mean, sd, and n
			vals_dict[capture_names[j]] = norm_signals

	# Mgic merging shit
	# Read in the merge groups and merge the items
	stats_dict = {}
	with open(replicate_groups, "r") as csvfile2:
		merge_lines_raw = csvfile2.readlines()
	merge_lines = []
	for line in merge_lines_raw:
		if len(line) > 1:
			merge_lines.append(line)
	print(merge_lines_raw)
	print(merge_lines)

	for line in merge_lines:
		merge_items = []
		items = line.split(",")
		for item in items:
			#print(item)
			for particle_signal in vals_dict[item.replace("\n", "")]:
				merge_items.append(particle_signal)
		stats_dict[items[0].replace("\n", "")] = {}
		stats_dict[items[0].replace("\n", "")]["mean"] = np.mean(merge_items)
		stats_dict[items[0].replace("\n", "")]["sd"] = np.std(merge_items)
		stats_dict[items[0].replace("\n", "")]["n"] = len(merge_items)

	# Write out
	g = open(no_ext(in_csv)+"_summary.csv", "w")
	g.write("Name,Mean Norm. Signal,Stdev,N\n")
	for capture in stats_dict:
		g.write(str(capture)+","+str(stats_dict[capture]["mean"])+","+str(stats_dict[capture]["sd"])+","+str(stats_dict[capture]["n"])+"\n")
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
	if len(sys.argv) == 3:
		main(str(sys.argv[1]), str(sys.argv[2]))
	else:
		print("Check argument inputs: foo.py inCSV inReplicateGroups\n\n")
		exit()