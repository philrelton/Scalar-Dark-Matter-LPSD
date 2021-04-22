# Code to run the LPSD algorithm

## Basic Documentation:
Most documentation for this code can be found here: https://gitlab.aei.uni-hannover.de/geoq/lpsd
Changes to this code include a basic checkpointing method to allow for longer analyses. This
does not affect run methods. However, result files are now saved to file periodically throughout
the calculation.

The primary change made was to add a new option allowing for simple parallelisation. This change 
is as follows:
- The original option `-n`, allowing for the definition of the total number of desired frequency
bins, has been replaced with a new option `-J`. This now sets the total number of frequencies bins 
to be calculated across the entire spectrum.
- The option `-n` now describes the number of frequencies to be calculated in each job. 
- The new option `-N` indicates which batch of frequencies to calculate.

We also made some modifications to enforce the calculation of the requested number of frequency
bins at the cost of reduced number of averages.

A final change was to slightly speed up reading of data files. The only possible format of file 
is now a two column text file with time in the first column and the data value in the second.

### Example:
In our analysis we wished to calculate bins between 50 and 8192 Hz at a resolution such that 
each bin had width 10^-6 of the central frequency. We therefore require 5098893 frequency bins.

To perform this analysis in a time efficient manner we apply these options as:
- `-J 5098893` to set the total number of frequency bins
- `-n 5000` to set the number of calculations in each batch
- `-N X` where X has value between 0 to 1019 for each job

This will result in 1020 jobs, the first 1019 producing an output file with 5000 lines of data, 
the final file will contain the final 3893 results.

NOTE: Due to the nature of the calculation, lower frequencies take longer to calculate. As such 
      these batches will take significantly longer than higher frequency batches.

### LPSD Algorithm originally described in:
Tröbs, M. and Heinzel, G., 2006. Improved spectrum estimation from digitized time series on a logarithmic frequency axis. Measurement, 39(2), pp.120-129.
https://www.sciencedirect.com/science/article/pii/S026322410500117X?casa_token=WcRBlyyEABYAAAAA:TqUNIcSN2qWlFMFr1eHROqnafKYUFQq14yDCpJX6S8PE593F9P5LSSOCL4AxL90fxb3PR9gHJw
