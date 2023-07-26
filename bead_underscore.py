import sys
import os
from os import listdir
from os.path import isfile, join
import shutil


def main():
	#find all tiff or tif files
	mypath = os.getcwd()
	onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
	tiff_files = []
	for f in onlyfiles:
		if (f[-5:] == ".tiff") or (f[-4:] == ".tif"):
			tiff_files.append(f)

	for tiff_file in tiff_files:
		shutil.move(tiff_file, no_ext(tiff_file).replace(" ", "_").replace("[", "_").replace("]", "_")+".tiff")


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