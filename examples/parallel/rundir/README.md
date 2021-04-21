# Example rundir for a parallelised LPSD calculation

## Run process
- Set run parameters in `submit.sh` script
    - These include the data file, the resolution the
frequency bounds and the number of frequencies per batch
- Submit via `sh submit.sh`
    - If you wish to resubmit a partially complete run then
`sh submit.sh` will do this
    - If you wish to completely restart a run then add the
`--restart` flag. WARNING: THIS WILL REMOVE ANY CURRENT RESULTS

The internall process of `submit.sh` is:
- All run variables are set
- Run variables are passed to the python script `Nfreqs.py`. This
calculates the number of frequency bins required at that resolution
- From this value, and the number of frequency bins desired in a batch,
it calculates the number of batches to run
- All values are then added to the `parallel.sh` script
- A condor dag is constructed to produce each batch
- The dag is then submitted to condor.
- All results appear in the `outdir/` directory
- Run files, including logs and errors are stored in the `dag_output/`
directory for each job.

This requires the `HTCondor` scheduler. Other schedulers will most likely
require adapations to this in order to run correctly.
