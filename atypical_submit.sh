#!/bin/bash

#snakemake --snakefile validate_atypicalwindows.smk \
#    --jobs 200 \
#    --rerun-incomplete \
#    --latency-wait 60 \
#    --cluster-config submit.json  \
#    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"


#!/bin/bash

for i in {1..3}
do
snakemake --snakefile validate_atypicalwindows.smk \
    --jobs 100 \
    all --batch all=${i}/3 \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"
done
