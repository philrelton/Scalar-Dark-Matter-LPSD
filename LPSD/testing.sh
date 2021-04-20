#!/usr/bin/env sh
#rm data_0.txt

tim=100
full=500
strt=4091.936
end=4195.524
n=500
#file=~/projects/geo-dark-matter/timeseries/testing_TS/1000s_TS_modified.txt
#out=modified.txt
file=../timeseries/testing_TS/timeseries_test.txt
out=modified.txt

#Declare iteration flag variable
for x in "$@"; do
    case $x in
        -i|--iteration)
        i="$2"
        shift 2
        ;;
    esac
done
i=0
./lpsd-exec \
	-A 2 \
	-b 0 \
	-e ${tim} \
	-f 16384 \
	-h 0 \
	-i ${file}\
	-l 20 \
	-n ${n} \
	-o ${out} \
	-r 0 \
	-s ${strt} \
	-t ${end} \
	-T \
	-w -2 \
	-p 238.13 \
	-x 1 \
	-N ${i} \
	-J ${full} \

#-i ~/projects/geo-dark-matter/testing/10000s_TS.txt \

# If you forget a variable (N or J in particular) everything will break and fill up the out files with 1000000000000 lines
