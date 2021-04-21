#!/usr/bin/env python3

import argparse
import numpy as np
from math import factorial

# Declare argparse variables
parser = argparse.ArgumentParser()
parser.add_argument("-m", "--fmin", default=50., dest="fmin", type=float)
parser.add_argument("-x", "--fmax", default=8192., dest="fmax", type=float)
parser.add_argument("-r", "--resolution", default=1e-6, dest="resolution", type=float)
parser.add_argument("-t", "--TSlength", default=5e4, dest="TSlength", type=float)

opts = parser.parse_args()

def lpsd_variables(fmin=50., fmax=8192., resolution=1e-6, TSlength=5e4):
    fsamp = 16384.
    Ndata = TSlength * fsamp
    overlap = 0.2
    
    g = np.log(fmax) - np.log(fmin)
    Nbins = int(1 + g/(np.log(fmin + fmin*resolution) - np.log(fmin)))
    jj = np.arange(Nbins, dtype=int)
    binwidth = fmin * np.exp(jj * g / float(Nbins - 1)) * (np.exp(g / float(Nbins - 1)) - 1)

    min_segment_length = np.around(float(fsamp) / binwidth).astype(int)[0]
    return Nbins, min_segment_length

Nbins, min_segment = lpsd_variables(fmin=opts.fmin, fmax=opts.fmax, resolution=opts.resolution, TSlength=opts.TSlength)

print(Nbins)

# Use these lines to check your time series is long enough
'''
if Ndata - min_segment <= 0:
    print(0)
    print("Your data is not long enough.")
    minTS = min_segment / 16384.
    print(f"You need a minimum timeseries of length: {minTS} seconds")
'''
