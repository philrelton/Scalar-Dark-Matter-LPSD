#!/usr/bin/env python3

import argparse
import numpy as np
from falsesignal import FalseSignal
import generalfunctions as GF

# Read in command line arguments
default_outdir =

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-o', '--outdir', default=default_outdir,
                    help="Output directory", type=str)
parser.add_argument('-m', '--modind', default=1,
                    help="Master Modification indicator", type=int)
opts = parser.parse_args()

# Check outdir ends correctly
if opts.outdir[-1] is not "/":
    opts.outdir += "/"

# Pick random number of signals: Integer
Nsignals = int(np.random.randint(4, 6, 1)[0])

#Create a False Signal object for each signal
signals = []
for idx in range(Nsignals): 
    signals.append(
		   FalseSignal(
		      frequency=float(np.random.uniform(200, 7000, 1)[0]),
		      amplitude=float(np.random.uniform(2e-22, 2e-21, 1)[0]),
		      phase_seed=int(np.random.randint(0, 2**32, 1)[0]),
		      Nfreqs=500,
		      FWHM=1e-6,
		      day=None,
		      leap_year=False
		      )
		   )
print("Signal parameters created")

# Save signal parameters: JSON file
GF.save_to_JSON(signals=signals,
             filename=f"{opts.outdir}MASTER_signals{opts.modind}")

