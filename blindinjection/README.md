# Job structure to create a time series with one or more injected signals

## `falsesignal.py`
Python script containing method to generate signal parameters for a dark
matter curve. Generates amplitdues along a defined curve shape with
frequencies and random phases assigned to each point on the curve.

## `generalfunctions.py`
Provides extra, general use, functions for signal creation

## `injection_dag.py`
Generates the submission files for condor

## `injection.sub`
Submits the jobs to condor
Currently set to generate the full dark matter curve. If a single sine
wave signal is required switch `executable` to `signal_generation_singlesine.py` 

## `parameter_generation.py`
Generates random parameters for a random number of signals
Creates a file in the outdir called `MASTER_signalsX.json`, where `X` is
an identifier for the injection. All injected time series will contain this
identifier in the file name

## `signal_generation.py`
Generates sine waves for each frequency across the curve across the full
period of the time series supplied.
Creates a new JSON containing the doppler shifted signal parameters
Creates a new time series containing the combination of the injected 
signals and the original time series

## `signal_generation_singlesine.py`
Generates a single sine wave at the desired frequency across the full
period of the time series
Creates a new JSON containing the doppler shifted signal parameters
Creates a new time series containing the combination of the injected 
signals and the original time series

## `submit.sh`
Submission for the main condor job. Injection indicator, out directory 
and input data sets are supplied here.
