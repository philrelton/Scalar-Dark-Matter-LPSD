#!/usr/bin/env python3

import numpy as np
import argparse
from tqdm import trange
import gwpy
from falsesignal import FalseSignal
import generalfunctions as GF  

# Read in command line arguments
default_outdir = 
default_TSloc = f"{default_outdir}GEO_TStest_s_1150300728_l_100.txt"

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', '--datafile', default=default_TSloc,
                    help="Original time series", type=str)
parser.add_argument('-o', '--outdir', default=default_outdir,
                    help="Output directory", type=str)
parser.add_argument('-m', '--modind', default=1,
                    help="Master modification indicator", type=int)
opts = parser.parse_args()

# Check outdir ends correctly
if opts.outdir[-1] is not "/":
    opts.outdir += "/"

# Load JSON of parameters: Dictionary
signals_dict = GF.load_from_JSON(
                   filename=f"{opts.outdir}/MASTER_signals{opts.modind}")

# Load real TS
print("Loading original timeseries")
original_TS = np.loadtxt(opts.datafile).T
print("Original time series loaded")

# Check real TS for day and leap year
gps_time = GF.name_split(opts.datafile)
MASTER_day, MASTER_leap_year = GF.day_year(gps_time)

# Create false signals from these: List of FS objects
signals = []
for signal in signals_dict.values(): 
    signals.append(
		   FalseSignal(
		      frequency=signal['frequency'],
		      amplitude=signal['amplitude'],
		      phase_seed=signal['phase_seed'],
		      Nfreqs=signal['Nfreqs'],
		      FWHM=signal['FWHM'],
		      day=MASTER_day,
		      leap_year=MASTER_leap_year
		      )
		   )
print("Time specific signals created")

# Save signals to JSON
filename = f"{opts.outdir}Time_specific_signals_{gps_time}_{opts.modind}"
GF.save_to_JSON(signals=signals, filename=filename)
print("Time specific signal parameters saved")

# Create signal
print("Adding signals to data...")
for signal in signals:
    for idx in trange(signal.Nfreqs):
        original_TS[1] = np.add(
	                original_TS[1],
                        signal["amplitudes"][idx]
                            * np.sin(
                                (2 * np.pi
                                 * signal["frequencies"][idx]
                                 * original_TS[0])
                                + signal["phases"][idx]
                                    )
                           )
print("False Signal(s) Created")

# Save new TS
print("Saving new time series...")
datafile_split = opts.datafile.rsplit("/", 1)[-1].split(".")[0]
outfile = f"{opts.outdir}{datafile_split}_{opts.modind}.txt"
np.savetxt(outfile, original_TS.T)
print("New timeseries saved.")
