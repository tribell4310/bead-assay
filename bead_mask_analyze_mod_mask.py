"""

bead_mask_analyze_mod_mask.py

Analysis scripts for 16-bit tiff images of fluorescent beads.
If using in multichannel analysis mode, hardcoded to look for a multichnnel.csv file that
lists the masking channel and analysis channel names, i.e.:
DAPI,FITC,Cy3

Adapting to have edge detection and internal circular masking.


"""

import sys
import imageio
import numpy as np
import matplotlib as mpl
mpl.use("Agg")
mpl.rc("figure", max_open_warning = 0)
import matplotlib.pyplot as plt
import math
import cv2 as cv
import skimage
from skimage.feature import peak_local_max
from scipy import ndimage as ndi
import os
import bead_mask_analyze_cython_modules as my_cy
from datetime import datetime
import statistics


def main(multichannel_flag, inFiles):
	start_time = datetime.now()

	# Create a logger to hold a dynamic copy of stdout
	log_stdout = open("run.out", "w")

	# Open the output file to write to.
	logfile = open(no_ext(inFiles[0])+"_log.txt", "w")
	logfile.write("Run started at "+str(start_time)+"\n\n")
	g = open(no_ext(inFiles[0])+"_E_analysis.csv", "w")
	g.write("Micrograph Name, Particle Number, Mask Size (px), Mean Signal, Median Signal, Min Signal, Max Signal, Std Dev\n")

	# Welcome message
	log_and_print("\nInput micrographs: "+str(len(inFiles))+"\n", log_stdout)

	# Multichannel check
	if multichannel_flag == 1:
		log_and_print("Performing analysis on multiple channels:", log_stdout)
		# Check for a file called multichannel.csv
		# If it exists, then parse out the mask_channel and analysis_channels.
		if os.path.isfile("multichannel.csv"):
			with open("multichannel.csv", "r") as f:
				contents = f.read().split(",")
				mask_channel = contents[0]
				analysis_channels = contents[1:]
			log_and_print("\tDetect beads on channel: "+mask_channel, log_stdout)
			log_and_print("\tApply masks to channels: "+ ", ".join(analysis_channels)+"\n", log_stdout)
		else:
			log_and_print("No multichannel.csv file was detected in this directory.  Exiting...", log_stdout)
			exit()
	else:
		log_and_print("Performing analysis only on the specified micrograph.\n", log_stdout)

	# Standard per-micrograph loop
	global_counter = 1
	problem_counter = 0
	global_error_flag = False
	for inFile in inFiles:
		# Load image
		log_and_print("\nLoading micrograph ("+str(global_counter)+"/"+str(len(inFiles))+"): "+inFile, log_stdout)
		image = imageio.imread(inFile)
		log_and_print("Image dimensions: "+str(image.shape[0])+" x "+str(image.shape[1])+" px ("+str(image.dtype)+")", log_stdout)
		#plt.imshow(image, cmap=plt.cm.gray)
		#plt.savefig(no_ext(inFile)+"_A_inputimage.png")
		#plt.clf()

		# Create a binary mask 
		image2, continue_flag = otsu_threshhold(image)
		if continue_flag == False:
			log_and_print("\n\tThis image doesn't appear to have any beads.  Moving on...\n", log_stdout)
			global_counter += 1
			continue

		# Generate a background mask and image
		bg_mask = background_mask_dt(image2)
		bg_image = my_cy.apply_mask(image, bg_mask)

		# Analyze the background image.
		filtered_bg_image, continue_flag_bg = filter_bg_image(bg_image, inFile, g)
		if continue_flag_bg == False:
			logfile.write("Micrograph skipped: "+inFile+"\n")
			global_error_flag = True
			global_counter += 1
			problem_counter += 1
			continue

		# Segmentation - resolve individual features.
		particles, particle_mask_sizes = simple_watershed(image2, image, inFile, bg_mask, bg_image, filtered_bg_image)
		if len(particles) == 0:
			log_and_print("\n\tNo beads passed size and edge filtering.  Moving on...\n", log_stdout)
			global_counter += 1
			continue

		# ##### Use the binary bask to create internal circular mask #####
		masked_particles, masked_particle_mask_sizes = apply_circular_mask(particles, particle_mask_sizes)
		if masked_particles == False:
			log_and_print("Problem with micrograph: "+inFile+"...skipped.\n", log_stdout)
			global_error_flag = True
			global_counter += 1
			problem_counter += 1
			continue

		# Multiply the original image by the binary mask to isolate and analyze individual particles.
		analyze_particle_images(inFile, g, image, masked_particles, masked_particle_mask_sizes)
		draw_originals(inFile, g, image, particles, particle_mask_sizes, masked_particles)

		global_counter += 1

		# If in multichannel mode, apply the same masks to a different image and report out
		if multichannel_flag == 1:
			for channel in analysis_channels:
				# Load the channel micrograph
				this_file = inFile.replace(mask_channel, channel.strip())
				log_and_print("Loading associated micrograph: "+str(this_file), log_stdout)
				channel_image = imageio.imread(this_file)
				
				# Apply the background mask to the new channel
				channel_bg_image = my_cy.apply_mask(channel_image, bg_mask)

				# Analyze the background image.
				filtered_bg_image, continue_flag_bg = filter_bg_image(channel_bg_image, this_file, g)
				if continue_flag_bg == False:
					logfile.write("IMAGE PARTICALLY PROCESSED: "+inFile+"\n")
					global_error_flag = True
					probem_counter += 1
					continue

				# Multiply the original image by the binary mask to isolate and analyze individual particles.
				analyze_particle_images(this_file, g, channel_image, masked_particles, masked_particle_mask_sizes)

	# Close the output file
	g.close()
	logfile.close()

	# Report back
	end_time = datetime.now()
	delta = end_time - start_time
	log_and_print("Micrographs processed:\t"+str(global_counter-1-problem_counter), log_stdout)
	log_and_print("Micrographs failed:\t"+str(problem_counter), log_stdout)
	log_and_print("Total run time:\t"+str(delta), log_stdout)
	if global_error_flag == True:
		log_and_print("\n"+(60*"*")+"\nNOTE: At least one micrograph was incompletely processed.  Check the logfile.\n"+(60*"*")+"\n", log_stdout)

	# Kill the stdout logger
	log_stdout.close()	


def worst_case_find_center_of_mask(inMask):
	""" This is a worst-case algorithm for finding the center of a binary mask.  Very slow. """
	total_x = 0
	total_y = 0
	num_px = 0
	
	for i in range(0, inMask.shape[0]):
		for j in range(0, inMask.shape[1]):
			if inMask[i][j] == 1:
				total_x += i
				total_y += j
				num_px += 1
	
	center_x = int(total_x / num_px)
	center_y = int(total_y / num_px)
	
	return (center_x, center_y), num_px


def apply_circular_mask(particles, particle_mask_sizes):
	masked_particles = []
	masked_particle_mask_sizes = []
	for particle in particles:
		# Apply distance transform and get maximum
		test_distance = 1
		total_peaks = 0
		problem_iters = 0
		dist_transform = ndi.distance_transform_edt(particle)
		# This block ensures we only have one center position, even if there were errors in segmentation
		while total_peaks != 1:
			escape_flag = False
			problem_iters += 1
			if problem_iters > 10: # risky loop with explicit defuse
				escape_flag = True
			local_max_identities = peak_local_max(dist_transform, min_distance=test_distance)
			total_peaks = len(local_max_identities)
			test_distance += 1
			if escape_flag == True:
				total_peaks = 1
		if escape_flag == True:
			center, size = worst_case_find_center_of_mask(particle)
		else:
			center = local_max_identities[0]

		# Find three points on mask edge by scanning from left, right, and top
		three_points = []
		# Left
		for i in range(0, particle.shape[0]):
			if particle[i][center[1]] == 1:
				three_points.append((i, center[1]))
				break
		# Right
		for i in range(1, particle.shape[0]+1):
			if particle[-i][center[1]] == 1:
				three_points.append((particle.shape[0]-i, center[1]))
				break
		# Top
		for i in range(0, particle.shape[1]):
			if particle[center[0]][i] == 1:
				three_points.append((center[0], i))
				break

		# Test to make sure this worked - not guaranteed to work for badly segmented particles where overall mask is not circular
		if len(three_points) != 3:
			# Keep previous rough center measurement and assign an arbitrary radius that does not impinge on the edges
			try: # this try-catch logic prevents weird shapes from breaking the program
				limited_radius = max([particle.shape[0] - center[0] - 1, particle.shape[1] - center[1] - 1], 50)
			except:
				return False, False

		else:
			# Calculate the circle that intersects these three points
			center, radius = circle_three_points(three_points[0], three_points[1], three_points[2])

			# Define a circle with same center but 0.8x radius (-> ~0.6x total area)
			limited_radius = int(radius * 0.8)

		# Use this circle to generate a binary mask
		new_mask = np.zeros((particle.shape[0], particle.shape[1]), dtype="uint16")
		rr, cc = skimage.draw.circle(center[0], center[1], limited_radius, (particle.shape[0], particle.shape[1]))
		new_mask[rr, cc] = 1
		
		# Append to return list
		masked_particles.append(new_mask)
		masked_particle_mask_sizes.append(cv.countNonZero(new_mask))

	return masked_particles, masked_particle_mask_sizes


def circle_three_points(A, B, C):
	# Define the slopes and intercepts of the perpendicular bisectors of vectors AB and BC
	m_perpAB, b_perpAB = perp_bisect(A, B)
	m_perpBC, b_perpBC = perp_bisect(B, C)

	# Determine the intersection of the perpendicular bisectors, this is the *origin*
	origin = get_intersection(m_perpAB, b_perpAB, m_perpBC, b_perpBC)
	if origin == False:
		return origin, False

	# Find the euclidian distance between the origin and A, this is the *radius* 
	radius = statistics.mean([get_euclidian_distance(A, origin), get_euclidian_distance(B, origin), get_euclidian_distance(C, origin)])

	return origin, radius


def perp_bisect(A, B):
	# If the y-vals are the same, we get meaningless values.  Have to add a pseudocount.
	if A[1] == B[1]:
		B = (B[0], (B[1]+0.0001))
	if A[0] == B[0]:
		B = ((B[0]+0.0001), B[1])

	# Define slope and intercept for vector AB
	m = (B[1] - A[1]) / (B[0] - A[0])
	b = A[1] - (m * A[0])

	# Slope of perpendicular is inverse reciprocal
	m_perp = -1 / m

	# Find the midpoint of AB
	midpoint_x = (A[0] + B[0]) / 2
	midpoint_y = (A[1] + B[1]) / 2

	# Find the intercept that fulfills the midpoint and the known slope
	b_perp = midpoint_y - (m_perp * midpoint_x)

	return m_perp, b_perp


def get_intersection(m1, b1, m2, b2):
	x_int = (b2 - b1) / (m1 - m2)
	y_int = (m1 * x_int) + b1
	try:
		return (int(round(x_int, 0)), int(round(y_int, 0)))
	except:
		return False


def get_euclidian_distance(A, B):
	y1 = A[1]
	y2 = B[1]
	x1 = A[0]
	x2 = B[0]
	return math.sqrt((y2-y1)**2 + (x2-x1)**2)


def log_and_print(inStr, logfile):
	print(inStr)
	logfile.write(inStr+"\n")
	return logfile


def analyze_particle_images(inFile, g, image, particles, mask_sizes):
	particle_counter = 1
	# Loop over every individual particle mask
	for i in range(len(particles)):
		if isinstance(particles[i], str) == False:
			# First, apply the mask - this can already be accelerated
			particle_image = my_cy.apply_mask(image, particles[i])

			fig, axes = plt.subplots(ncols=1, figsize=(8, 8), sharex=True, sharey=True)
			#ax = axes.ravel()
			plt.imshow(particle_image, cmap=plt.cm.gray)
			plt.savefig(no_ext(inFile)+"_D_"+pad_to_int(particle_counter, 3)+"_particleimage.png")
			plt.clf()

			# Analysis
			try:
				num_white = mask_sizes[i]
				signals = []
				for i in range(0, particle_image.shape[0]):
					for j in range(0, particle_image.shape[1]):
						if particle_image[i, j] != 0:
							signals.append(particle_image[i, j])

				min_intensity = min(signals)
				max_intensity = max(signals)
				mean_intensity = np.mean(signals)
				median_intensity = np.median(signals)
				stdev_intensity = np.std(signals)
			except:
				break

			g.write(str(inFile)+", "+str(particle_counter)+", "+str(num_white)+", "+str(mean_intensity)+", "+str(median_intensity)+", "+str(min_intensity)+", "+str(max_intensity)+", "+str(stdev_intensity)+"\n")
		particle_counter += 1


def draw_originals(inFile, g, image, particles, particle_mask_sizes, masked_particles):
	particle_counter = 1
	position_counter = 0

	# Loop over every individual particle mask and write out original mask per_beadimages
	for i in range(len(masked_particles)):
		# Check whether particle was filtered
		if isinstance(masked_particles[i], str) == False:
			# First, apply the mask - this can already be accelerated
			particle_image = my_cy.apply_mask(image, particles[position_counter])

			# Write out image to file
			fig, axes = plt.subplots(ncols=1, figsize=(8, 8), sharex=True, sharey=True)
			#ax = axes.ravel()
			plt.imshow(particle_image, cmap=plt.cm.gray)
			plt.savefig(no_ext(inFile)+"_C_"+pad_to_int(particle_counter, 3)+"_particleimage_fullBeadSegment.png")
			plt.clf()

			position_counter += 1

		particle_counter += 1


def filter_bg_image(bg_image, inFile, g):
	# First, plot the pixel intensities across the non-zero pixels
	non_zero_pixels = 0
	non_zero_values = []
	for x in range(0, bg_image.shape[0]):
		for y in range(0, bg_image.shape[1]):
			if bg_image[x][y] != 0:
				non_zero_values.append(bg_image[x][y])
				non_zero_pixels += 1
	
	# Remove any two-sigma outliers
	filtered_values = []
	filtered_values_normalized = []
	filtered_pixels = 0
	bg_mean = np.mean(non_zero_values)
	bg_stdev = np.std(non_zero_values)
	lower_bound = bg_mean - (2*bg_stdev)
	upper_bound = bg_mean + (2*bg_stdev)

	for non_zero_value in non_zero_values:
		if (non_zero_value > lower_bound) and (non_zero_value < upper_bound):
			filtered_values.append(non_zero_value)
			filtered_values_normalized.append(non_zero_value)
			filtered_pixels += 1

	try:
		# Generate a version of the bg_image after filtration
		filtered_bg_image = my_cy.apply_filter_to_bg_image(bg_image, lower_bound, upper_bound)
		# Write out to reporting doc
		norm = 65535
		g.write(str(inFile)+", BG, "+str(filtered_pixels)+", "+str(np.mean(filtered_values)/norm)+", "+str(np.median(filtered_values)/norm)+", "+str(min(filtered_values)/norm)+", "+str(max(filtered_values)/norm)+", "+str(np.std(filtered_values)/norm)+", Comments: "+str(non_zero_pixels - filtered_pixels)+" pixels with two-sigma intensity below "+str(bg_mean - (2*bg_stdev))+" or above "+str(bg_mean + (2*bg_stdev))+" were removed.\n")
		continue_flag = True
	except:
		filtered_bg_image = False
		continue_flag = False

	# Remove micrograph from analysis if the background has less than 2000 non-zero pixels (overcrowded image)
	if filtered_pixels < 2000:
		continue_flag = False
	
	return filtered_bg_image, continue_flag


def background_mask_dt(inMask):
	# Define some static variables to avoid unnecessary calls
	mask_shape_x = inMask.shape[0]
	mask_shape_y = inMask.shape[1]

	# Create the new mask object
	bg_mask_negative = np.copy(inMask) # white circles on black bg
	bg_mask = np.invert(bg_mask_negative) # black circles on white bg

	# Calculate a distance transform of bg_mask
	bg_mask_dt = ndi.distance_transform_edt(bg_mask)

	# Pass any pixels within euclidian distance 70 of a black pixel to the bit flipper
	# or anything within 70px of the edge
	for x in range(0, mask_shape_x):
		for y in range(0, mask_shape_y):
			flip_flag = False
			value = bg_mask_dt[x][y]
			edge_border = 150
			if value <= 140: # use for agarose resin
			#if value <= 70: # use for sepharose resin
				flip_flag = True
			if (x < edge_border) or (x > (mask_shape_x - edge_border)):
				flip_flag = True
			if (y < edge_border) or (y > (mask_shape_y - edge_border)):
				flip_flag = True

			if flip_flag == True:
				flip_px_to_zero(bg_mask, x, y)

	return bg_mask


def flip_px_to_zero(mask, x, y):
	mask[x, y] = 0
	return mask


def otsu_threshhold(inImg):
	# Actually perform threshholding.
	otsu_threshold, image_result = cv.threshold(inImg, 0, 65535, cv.THRESH_BINARY + cv.THRESH_OTSU)
	
	# If the threshhold is close to the max intensity, this image likely doesn't have any beads... discard.
	if float(otsu_threshold) / float(inImg.max()) > 0.62:
		continue_flag = False
	else:
		continue_flag = True

	return(image_result, continue_flag)


def simple_watershed(img, orig_img, inFile, bg_mask, bg_image, filtered_bg_image):
	# Calculate distance transform of binary mask and find local maxima
	dist_transform = ndi.distance_transform_edt(img)
	local_max_identities = peak_local_max(dist_transform, min_distance=1)
	local_max_boolean = np.zeros((img.shape[0], img.shape[1]), dtype="uint16")
	for local_max in local_max_identities:
		local_max_boolean[local_max[0]][local_max[1]] = 1

	# Segment each local maximum into a marked object and extract its properties
	markers, _ = ndi.label(local_max_boolean)
	segmented = skimage.segmentation.watershed(65535-dist_transform, markers, mask=img, compactness=0.005)
	object_labels = skimage.measure.label(segmented)
	properties = skimage.measure.regionprops(object_labels)

	# Output plot
	fig, axes = plt.subplots(ncols=7, figsize=(21, 3), sharex=True, sharey=True)
	ax = axes.ravel()
	ax[0].imshow(orig_img, cmap=plt.cm.gray)
	ax[0].set_title('Input Image')
	ax[1].imshow(img, cmap=plt.cm.gray)
	ax[1].set_title('Binary Mask')
	ax[2].imshow(-dist_transform, cmap=plt.cm.gray)
	ax[2].set_title('Distance transform')
	ax[3].imshow(segmented, cmap=plt.cm.nipy_spectral)
	ax[3].set_title('Separated objects')
	ax[4].imshow(bg_mask, cmap=plt.cm.gray)
	ax[4].set_title('Background mask')
	ax[5].imshow(bg_image, cmap=plt.cm.nipy_spectral)
	ax[5].set_title('Background image')
	ax[6].imshow(filtered_bg_image, cmap=plt.cm.nipy_spectral)
	ax[6].set_title('Filtered bg image')
	for a in ax:
	    a.set_axis_off()
	fig.tight_layout()
	plt.savefig(no_ext(inFile)+"_B_binarymask.png")
	plt.clf()

	# Start by filtering anything truly tiny - less than 2500 px2 for agarose resins
	val_properties = []
	for i in range(0, len(properties)):
		if properties[i].area > 2500: #use for agarose resins
		#if properties[i].area > 200: #use for sepharose resins
			val_properties.append(properties[i])
	
	# Merge anything where the centroids are too close.
	merge_container = []
	for i in range(0, len(val_properties)):		
		for j in range(i+1, len(val_properties)):
			if euclidian_distance(val_properties[i].centroid, val_properties[j].centroid) < 100:
				merge_container.append((i, j))

	merge_containter_items = []
	for item in merge_container:
		merge_containter_items.append(item[0])
		merge_containter_items.append(item[1])

	# Create index list groups.
	merge_groups = []
	merge_group_contents = []
	for i in range(0, len(val_properties)):
		# Update list merge_group_contents with the current contents of merge_groups
		for item in merge_groups:
			for thing in item:
				if thing not in merge_group_contents:
					merge_group_contents.append(thing)

		if i not in merge_containter_items:
			merge_groups.append([i])

		else:
			if i not in merge_group_contents:
				friends = [i]
				# Find all its pairs
				for j in range(0, len(merge_container)):
					if (merge_container[j][0] == i):
						friends.append(merge_container[j][1])
					elif (merge_container[j][1] == i):
						friends.append(merge_container[j][0])
				merge_groups.append(friends)

	# Actual implementation of block merging
	# Pass all relevant pixels from each item in merge_groups to the binary mask creator
	
	val_mask_sizes = []
	for i in range(len(val_properties)):
		val_mask_sizes.append(val_properties[i].area)

	outputs = []
	outputs_mask_sizes = []
	counter = 0
	for merge_group in merge_groups:
		counter += 1
		merge_area = 0
		# Make the new array
		particle_mask = np.zeros((img.shape[0], img.shape[1]), dtype="uint16")
		for properties_index in merge_group:
			for pixel in val_properties[properties_index].coords:
				particle_mask[(pixel[0], pixel[1])] = 1
			merge_area += val_properties[properties_index].area
		outputs.append(particle_mask)
		outputs_mask_sizes.append(merge_area)

	# Here need to add a final edge detection check, and only pass outputs that pass that check.
	# Can't be within some px of an edge
	pass_outputs = []
	pass_outputs_mask_sizes = []
	for i in range(0, len(outputs)):
		if test_edge(outputs[i], 20) == True:
			pass_outputs.append(outputs[i])
			pass_outputs_mask_sizes.append(outputs_mask_sizes[i])
			
	return pass_outputs, pass_outputs_mask_sizes


def test_edge(inMask, thresh):
	pass_marker = True
	square_length = inMask.shape[0]
	# Scan edges - this assumes that the image is square.
	for i in range(0, square_length):
		for j in range(0, thresh):
			if (inMask[j][i] == 1) or (inMask[square_length-j-1][i] == 1) or (inMask[i][j] == 1) or (inMask[i][square_length-j-1] == 1):
				pass_marker = False 
				break

	return pass_marker


def euclidian_distance(pt1, pt2):
	return math.sqrt(((float(pt1[1])-float(pt2[1]))**2) + ((float(pt1[0])-float(pt2[0]))**2))


def pad_to_int(inInt, n):
	inStr = str(inInt)
	while len(inStr) < n:
		inStr = "0"+inStr
	return inStr


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
		main(int(sys.argv[1]), sys.argv[2:])
	else:
		print("Check argument inputs: foo.py multichannel inFile1 inFile2 inFile3 ...\n\n for multichannel, 0=just process the listed micrograph, 1=use multichannel.csv to apply these masks to other images")
		exit()
