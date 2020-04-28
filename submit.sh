#!/bin/bash

for i in {1..4}; do
snakemake --jobs 200 \
    --batch all=${i}/4 \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"; done

#snakemake --jobs 200 \
#    --batch all=2/4 \
#    --rerun-incomplete \
#    --latency-wait 60 \
#    --cluster-config submit.json  \
#    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"

#snakemake --jobs 200 \
#    --batch all=3/4 \
#    --rerun-incomplete \
#    --latency-wait 60 \
#    --cluster-config submit.json  \
#    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"

#snakemake --jobs 200 \
#    --batch all=4/4 \
#    --rerun-incomplete \
#    --latency-wait 60 \
#    --cluster-config submit.json  \
#    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"

