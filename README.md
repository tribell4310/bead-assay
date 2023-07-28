# bead-assay

A quantitative, hypersensitive, and lipid-compatible IP assay.

![Capture_1_-_Position_1__81__DAPI_26_B_binarymask](https://github.com/tribell4310/bead-assay/assets/67428134/7103ba8a-5de7-49fe-b770-34d6ff0e1a79)

# Approach

The bead assay is a twist on a traditional immunourification (IP) assay.  In a traditional IP, one protein (the "bait") is bound to a resin and the amount of a second protein (the "prey") is measured, usually by SDS-PAGE or immunoblotting.  Though robust, this method is not particularly quantitative, and is limited to looking at protein preys.

We extend the IP assay format by processing bound prey signal using fluorescence microscopy and computational analysis.  This allows us to measure the binding of any type of prey (proteins, ligands, glycans, lipids, etc.) so long as it is fluorescently tagged.  We can also exploit the fact that individual agarose or sepharose resin beads are spherical to eliminate air bubbles, cellular debris, or other artifacts that could interfere with sensitive measurements.

![image](https://github.com/tribell4310/bead-assay/assets/67428134/2bd66b5f-925a-4e2a-ab61-5114d9580863)

This approach will soon be published in our new preprint article on biorxiv!  In the meantime, please inquire directly for detailed protocols: bell@molbio.mgh.harvard.edu.

## The IP

This protocol is amenable to any existing IP protocol so long as the prey molecule is fluorescently tagged.  We regularly measure BODIPY-cholesterol binding to DDM-solubilized Strep-tagged Prominin-1 using an anti-Strep resin.  

We also trace label our resin with mTagBFP2 prior to binding the bait molecule.  We have a Strep- and His- tagged construct and a FLAG and HA tagged construct for bacterial expression that will soon be available at Addgene.  In the meantime, please feel free to inquire directly.  This trace label is used to identify the resin beads and to normalize fluorescence signals during computational processing.

Once your protein is bound to the resin is bound and washed (just prior to where you would normally elute), you are ready to begin the customized detection pipeline.

## Data collection

 1. Resuspend your washed resin in 40-50 uL of wash buffer.
 2. Using a cut pipette tip, transfer the beads to a glass microscope slide and gently cover with a piece of cover glass.  Try to avoid air bubbles becoming trapped between the slide and cover glass.
 3. Collect a montaged series of 80-120 images across the cover glass.  Image on the blue (DAPI) channel to detect beads and on a channel compatible with your fluorophore of interest.
 4. Export your images as 16-bit TIF images, using the following nomenclature.  To be compatible with the processing scripts, it is only strictly required that the channel name appear in the name directly before the file extension.

`experimentName_montagePositionNumber_channelName.tiff`

Some microscope software produces default filenames with unix-disallowed characters, or a `.tif` instead of a `.tiff`	file extension.  For convenience, the `bead_underscore.py` script is included that corrects these issues.

## Data processing

The data processing scripts require python 3.5 or higher, preferably installed in a dedicated virtual environment.  The scripts should be placed into a dedicated working directory. Package dependencies are:
 - cython
 - csv
 - math
 - imageio
 - numpy
 - matplotlib
 - cv2
 - skimage
 - scipy
 - datetime
 - statistics

In addition, your system will need:
 - C compiler (gcc on Linux systems or equivalent for MacOS).
 - imagemagick

### Pre-compile custom functions

These scripts accelerate image analysis by pre-compiling certain functions.  To compile, navigate to the working directory and run the following command:

`python build.py build_ext --inplace`

After compiling, a new file with a `.so` extension should appear.  This is your pre-compiled set of analysis functions.

### Split multi-page tiff files into single-frame tiffs

In our microscope setup, out exported tiff files are multi-page, meaning that each position of the slide we acquired on is a frame in a single tiff file.  We first need to split this multi-page tiff file into individual files.  You can do this with the provided `magick_convert.sh` convenience script.  On Linux, this would look like:

`bash magick_convert.sh`

This script will generate a single-page tiff for every frame in your multi-page tiffs.  Before proceeding, remove the multi-page tiffs from your working directory.

### Configure parameter files
The analysis scripts require two parameter files, alled `multichannel.csv` and `replicates.csv`.  The `multichannel.csv` file needs to contain the channels you want to measure fluorescence on, with the first channel being the one used for bead detection.  Because we image our BFP bead marker on the DAPI channel and image bound BODIPY-cholesterol on the FITC channel.  Thus, our multichannel file would look like:
> DAPI,FITC

The `replicates.csv` file defines constant portions of the names of the image files to look at.  For example, if some tiff image files look like this:
> Capture_1_DAPI_1<br>
> Capture_1_FITC_1<br>
> Capture_1_DAPI_2<br>
> Capture_1_FITC_2<br>
> Capture_2_DAPI_1<br>
> Capture_2_FITC_1<br>
> Capture_2_DAPI_2<br>
> Capture_2_FITC_2

...then the replicates file should look like:
> Capture_1<br>
> Capture_2

### Analyze and segment images of beads

Run the following command:

`python bead_analyze_everything_batches.py numProc channel multichannelFlag`

|Parameter|Type|Description|Example from data above|
|--|--|--|--|
| numProc | int | Number of parallel threads to use for processing (each thread will use up to four CPU cores at a time). | Anywhere from 1-5 depending on system resources. |
| channel | string | Fluorescence channel to detect beads on. | `DAPI` |
| multichannelFlag | boolean (0, 1) | Whether to quantify fluorescence on multiple channels (`0` only quantifies the detection channel, `1` quantified the detection channel and the other channels in `multichannel.csv`). | `1` |


This script will create a file called RUNME.sh.  Execute this file by running the command:

`bash RUNME.sh`

The script will process the images, and each thread will print a message to the terminal when it terminates successfully.  Keep an eye on your system resources to know when all threads have finished.

### Clean-up from automatic processing

Concatenate the outputs from multiple threads using the command:
`python bead_concat_csv.py`

This will produce a file called `collate_analysis.csv` containing all of the quantification data.

Next, convert the analysis spreadsheet into a human-readable form using the command:
`python bead_analyze_csv_output_mask.py collate_analysis.csv [channels]`

|Parameter|Type|Description|Example from data above|
|--|--|--|--|
| [channels] | stings, separated by spaces | Channels that quantification was performed on.  These are identical to the channels in you `multichannel.csv` file. | `DAPI FITC` |

Finally, run the following command to package the script outputs into a usable format.
`bash bead_cleanup.sh`

This will produce a directory called `copy_out` containing all of the relevant output data.  If working on a remote server, use sftp to copy these files to your local device.

## Manual Curation

To achieve the most sensitive possible measurements, it is necessary to manually review the images of resin beads and eliminate any that were incorrectly identified.  Images of beads are stored within the `copy_out/full_bead_segment_pngs` directory.  Real resin beads should look perfectly circular:

![image](https://github.com/tribell4310/bead-assay/assets/67428134/2aebc8c6-8365-413a-83d3-43055a3c0d74)


![image](https://github.com/tribell4310/bead-assay/assets/67428134/938439e0-f879-4757-a4d2-d4a4e8785424)


![image](https://github.com/tribell4310/bead-assay/assets/67428134/bc286df9-b7e2-4f54-b647-5034b4f67f91)


Artifacts are easy to spot by eye.

![image](https://github.com/tribell4310/bead-assay/assets/67428134/a2875f1f-1ddb-4b72-b3d9-2cf3aa69d08e)


![image](https://github.com/tribell4310/bead-assay/assets/67428134/098d7578-807f-468f-8c39-fde35d527f86)


Take note of the particle IDs for any objects that should be excluded from analysis.

Create a copy of the `collate_analysis_analyzed.csv` file.  I usually call mine `collate_analysis_analyzed_manCurate.csv` to indicate that I've curated it manually.  Open your new copy in any spreadsheet application and delete the rows corresponding to the artifacts that should be excluded from analysis.  Save your modified file as a csv and copy it back to the working directory.

## Final analysis

Finally, calculate summary statistics (averages, standard deviation, and replicate count) across the different beads for each of your conditions.  We provide a convenience script that reports out our preferred metric: the signal on each channel normalized against the bead detection channel.  In the examples above, the output values are background-corrected FITC/DAPI ratio.

`python bead_calc summary_stats.py collate_analysis_analyzed_manCurate.csv replicates.csv`

This produces a final output called `collate_analysis_analyzed_manCurate_summary.csv` containing your final collated data.

If you'd like to use a different measurement, all of the per-bead data (with and without background correction) is retained in the `collate_analysis_analyzed_manCurate.csv` spreadsheet file.


## Questions

If you run into issues with this pipeline, please feel free to open an Issue in GitHub, and I'll get back to you as quickly as I can.  Thanks for your interest!


> Written with [StackEdit](https://stackedit.io/).
