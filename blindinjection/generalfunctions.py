#!/usr/bin/env python3

import json
import numpy as np
import gwpy

def save_to_JSON(signals, filename="TEST"):
    '''
    Saves signal parameters to JSON format
    '''
    signal_storage = {}
    for N, signal in enumerate(signals):
        parameters = {}
        parameters["frequency"] = signal.frequency
        parameters["amplitude"] = signal.amplitude
        parameters["phase_seed"] = signal.phase_seed
        parameters["Nfreqs"] = signal.Nfreqs
        parameters["FWHM"] = signal.FWHM
        parameters["day"] = signal.day
        parameters["leap_year"] = signal.leap_year
        signal_storage[f"signal{N}"] = parameters

    out_file = open(f"{filename}.json",'w+')
    json.dump(signal_storage, out_file, indent=4)
    out_file.close()
    print("Signal Parameters Saved")

def load_from_JSON(filename):
    '''
    Loads signal parameters
    '''
    with open(f"{filename}.json", "r") as json_file:
        signals_dict = json.load(json_file)
    json_file.close()
    return signals_dict

def name_split(string_name):
    '''
    Splits the file name to acquire the gps time
    '''
    return int(string_name.rsplit("/", 1)[-1].split("_")[-3])

def day_year(gps_time):
    '''
    Calculates the day of the year to apply seasonal
    doppler shift
    '''
    MASTER_day = int(gwpy.time.from_gps(gps_time).strftime('%j'))

    leap_years = [idx for idx in range(2000, 2024, 4)]

    if gwpy.time.from_gps(gps_time).year in leap_years:
        MASTER_leap_year = True
    else:
        MASTER_leap_year = False
    return MASTER_day, MASTER_leap_year
