"""
concat_csv.py

When running simultaneous threads, I get multiple score outputs that need to be put back together.
 - Finds csv files
 - Concatenates them
 - Outputs them to a new "collate_analysis.csv" file.

"""


import sys
import os


def main():
	cwd = os.getcwd()
	counter = 0
	names = []
	for root, dirs, files in os.walk("."):
		for filename in files:
			if filename[-4:] == ".csv":
				if "multichannel" not in filename:
					if "collate_analysis" not in filename:
						if "replicates" not in filename:
							if os.path.isfile(os.path.join(cwd, filename)):
								names.append(filename)
								counter += 1
	#print(names)

	print("\nDetected "+str(counter)+" csv files.")

	g = open("collate_analysis.csv", "w")
	# Check to see if collate_analysis is in names - if so remove it
	first_Flag = False
	for name in names:
		if name != "collate_analysis.csv":
			#print(name)
			f = open(name, "r")
			lines = f.readlines()
			#print(lines)
			if first_Flag == False:
				g.write(lines[0])
				first_Flag = True
			for i in range(1, len(lines)):
				g.write(lines[i])
			f.close()
	g.close()

	print("...done!\n")


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
	if len(sys.argv) == 1:
		main()
	else:
		print("Check argument inputs: foo.py")
		exit()
