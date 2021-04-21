import argparse
import os

# Default values
default_outdir =
default_datafile =

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-f', '--datafiles', default=default_datafile,
                    help="Data file locations", type=str)
parser.add_argument('-o', '--outdir', default=default_outdir,
                    help="Output directory", type=str)
parser.add_argument('-m', '--modind', default=1,
                    help="Modificaiton indicator", type=int)
opts = parser.parse_args()
opts.datafiles = opts.datafiles.split(",")

# Writes jobs to condor dag script
with open("injection.dag", "w") as f:
    for idx in range(len(opts.datafiles)):
        f.writelines([f"JOB INJECTION{idx} injection.sub\n"])
    f.writelines(["\n"])
    for idx, datafile in enumerate(opts.datafiles):
        f.writelines([f'VARS INJECTION{idx}' \
                      f' datafile="{datafile}"' \
                      f' outdir="{opts.outdir}"' \
                      f' modind="{opts.modind}"' \
                      f' jobind="{idx}" \n'])
f.close()

# Create job directories
if not os.path.isdir("dag_output/"):
    os.mkdir("dag_output/")

for idx in range(len(opts.datafiles)):
    os.mkdir(f"dag_output/file{idx}_mod{opts.modind}/")
