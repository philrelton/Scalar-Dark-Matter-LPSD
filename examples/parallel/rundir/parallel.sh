#!/usr/bin/env sh

#Declare iteration flag variable
for x in "$@"; do
    case $x in
        -i|--iteration)
        iter="$2"
        shift 2
        ;;
    esac
done

# Declare data file
# Must be two column text file: Time Data
filename=/path/to/file/filename.txt

# Declare number of seconds data in file
TSlength=50000

# Declare start and end frequencies
# These cover the width of the entire job not the sub-jobs
fmin=50
fmax=8192

# Declare 'Jdes': The total number of required frequency bins
# This can also be declared by argument
Nfreqs_full=5098893

# Declare the number of frequencies in each batch
# If a single job is required:
# Set i=0 and n=${full}
batch_size=5000

../../../LPSD/lpsd-exec \
	-A 2 \
	-b 0 \
	-e ${TSlength} \
	-f 16384 \
	-h 0 \
	-i ${filename} \
	-l 20 \
	-n ${batch_size} \
	-o outputs/data_${iter}.txt \
	-r 0 \
	-s ${fmin} \
	-t ${fmax} \
	-T \
	-w -2 \
	-p 238.13 \
	-x 1 \
	-N ${iter} \
	-J ${Nfreqs_full} \

# If you forget a variable (N or J in particular) everything will break and fill up the out files with 1000000000000 lines
