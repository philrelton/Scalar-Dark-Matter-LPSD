import numpy as np
from tqdm import trange
import os.path

run = 0
bash = f"ls run{run}/outputs/ | wc -l"
nfile = int(os.popen(bash).read())

full_freq = np.loadtxt(f"run{run}/outputs/data_0.txt").T[0]
full_psd = np.loadtxt(f"run{run}/outputs/data_0.txt").T[1]
full_ps = np.loadtxt(f"run{run}/outputs/data_0.txt").T[2]
av = np.loadtxt(f"run{run}/outputs/data_0.txt").T[3]

if nfile > 1:
    for i in trange(1, nfile):
        file_name = f"run{run}/outputs/data_{i}.txt"
        if os.path.isfile(file_name):
            temp = np.loadtxt(file_name)
            full_freq = np.append(full_freq, temp.T[0])
            full_psd = np.append(full_psd, temp.T[1])
            full_ps = np.append(full_ps, temp.T[2])
            av = np.append(av, temp.T[3])

concat = np.vstack((full_freq, full_psd, full_ps, av)).T
np.savetxt(f"outputs/power_spectrum_TS{run}.txt", concat, fmt="%1.6e")

