# Example run structure for a parallelised LPSD calculation

## Running a parallel job:
- Create a copy of `rundir/` called `run0/`
- Move to that run directory, modify `submit.sh` as explained in
the README.md file in `rundir/`
- Submit the job using: `sh submit.sh`
- Once the results are finished run `python3 combine_dag.py`
    - This will combine the individual result files into a single
      spectrum.
    - The run directory will need specifying in `combine_dag.py`
      i.e. `run1/`
