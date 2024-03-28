#!/bin/bash
#SBATCH -t 150:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE
#SBATCH --partition uri-cpu

module load uri/main
module load snakemake/6.10.0-foss-2021b

snakemake --unlock
snakemake --cluster "sbatch -t {cluster.time} -N {cluster.nodes} --ntasks-per-node {cluster.ntasks-per-node} --mem {cluster.mem} --partition {cluster.partition}" --latency-wait 120 --cluster-config config.yml -j 10

