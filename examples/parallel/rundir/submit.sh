# Check for restart flag
for x in "$@"; do
    case $x in
        -r|--restart)
        rest=1
        shift 2
        ;;
    esac
done

# Set and calculate run parameters
fmin=50
fmax=8192
resolution=1e-6
TSlength=50000
filename=/path/to/file/filename.txt
# Declare the size of each batch
batch_size=5000

# Calculate the required number of bins for this resolution
Nfreqs_full=`python3 Nfreqs.py --fmin ${fmin} --fmax ${fmax} --resolution ${resolution} --TSlength ${TSlength}`
# Be careful here. If your time series is not long enough then the low frequencies will fail at high resolution
# There are a few lines to check in Nfreqs.py if you are unsure

# Calculate the number of required batches
# The final batch will contain between one and batch_size bins
Nbatches=`python -c "import numpy as np; print(np.ceil(${Nfreqs_full} / ${batch_size}).astype(int))"`

# Replace run parameters in parallel.sh file
sed -i "s/fmin=.*/fmin=${fmin}/" parallel.sh
sed -i "s/fmax=.*/fmax=${fmax}/" parallel.sh
sed -i "s/TSlength=.*/TSlength=${TSlength}/" parallel.sh
sed -i "s;filename=.*;filename=${filename};" parallel.sh
sed -i "s/Nfreqs_full=.*/Nfreqs_full=${Nfreqs_full}/" parallel.sh
sed -i "s/batch_size=.*/batch_size=${batch_size}/" parallel.sh

# Check if restart is desired and generate run files/dir
if [ -z "${rest}" ]; then
    if [ -f "lpsd.dag" ]; then
        echo "Continuing run"
    else
        echo "Starting new run"    
        mkdir outputs/
        python3 ../generate_dag.py --number ${Nbatches}
    fi
else
    echo "Restarting run"
    rm lpsd.dag*
    rm -r dag_output/*
    rm outputs/*
    python3 ../generate_dag.py --number ${Nbatches}
fi

# Submit job
condor_submit_dag lpsd.dag
