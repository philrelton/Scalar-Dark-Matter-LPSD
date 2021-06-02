#!/usr/bin/env sh

# Declare data file
# Must be two column text file: Time Data
filename=path/to/data/file.txt

# Declare number of seconds data in file
TSlength=100

# Declare start and end frequencies
fmin=50
fmax=8192

# Declare 'Jdes': The total number of required frequency bins
Nfreqs_full=1000

# Declare the number of frequencies in each batch
# If a single job is required:
# Set iter=0 and n=${full}
iter=0
batch_size=${Nfreqs_full}

# Declare the outup file
outfilename=result.txt

../../LPSD/lpsd-exec \
	-A 2 \
	-b 0 \
	-e ${TSlength} \
	-f 16384 \
	-h 0 \
	-i ${filename}\
	-l 20 \
	-n ${batch_size} \
	-o ${outfilename} \
	-r 0 \
	-s ${fmin} \
	-t ${fmax} \
	-T \
	-w -2 \
	-p 238.13 \
	-x 1 \
	-N ${iter} \
	-J ${Nfreqs_full} \

# If an argument is missed then the original code will request commandline input
# Be careful this does not occur in a job submitted to a compute node as this
# will quickly create large files with repeated requests.
