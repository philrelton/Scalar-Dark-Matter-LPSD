import argparse
import os
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('-n', '--number', default=100,
                    help="Number of jobs to generate", type=int)
opts = parser.parse_args()

with open("lpsd.dag", "w") as f:
    for idx in range(opts.number):
        f.writelines([f"JOB PSD{idx} parallel.sub\n"])
    f.writelines(["\n"])
    for idx in range(opts.number):
        f.writelines([f"VARS PSD{idx} iteration=\"{idx}\"\n"])

if not os.path.isdir("dag_output/"):
    os.mkdir("dag_output/")

for idx in range(opts.number):
    os.mkdir(f"dag_output/run_{idx}/")

