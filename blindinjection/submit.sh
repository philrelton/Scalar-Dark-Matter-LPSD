## Initialise parameters
# Set directory for outputs
outdir=
mkdir ${outdir}
# Set indicator for which set of injections this is
modind=1
# Set location of data files
# Data files must have the form AA_<GPSTIME>_<Period>.txt
# Data files must be text files with two columns: Time  Data
datadir=
file1=${datadir}GEO_TS0_s_1150300728_l_69000.txt
file2=${datadir}GEO_TS22_s_1167087798_l_69000.txt
datafiles=${file1},${file2}

# Create injection parameters
python parameter_generation.py \
    -o ${outdir} \
    -m ${modind}

# Make run files
python injection_dag.py \
    -f ${datafiles} \
    -o ${outdir} \
    -m ${modind}

# Submit jobs
condor_submit_dag injection.dag

