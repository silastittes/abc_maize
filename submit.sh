#!/bin/bash

snakemake --jobs 200 \
    --batch all=1/4 \
    --rerun-incomplete \
    --latency-wait 60 \
    --cluster-config submit.json  \
    --cluster "sbatch -p {cluster.p} --mem {cluster.mem} --cpus-per-task {cluster.cpus-per-task} --time {cluster.time} --job-name {cluster.name} -e {cluster.e} -o {cluster.o}"

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

