#!/bin/bash

for i in {1..1}
do
snakemake --jobs 200 \
    all --batch all=${i}/5 \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"
done

