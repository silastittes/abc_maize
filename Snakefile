import numpy as np
import random

#fixed params
loci = 6000000 ! HOW DOES THIS LOOK FOR REAL DATA??
sample = 26
mu = 3e-7
rr = 1.6e-6

#priors
n_draws = 100000

#random.seed(214125) #random seed per sim is reported, so no reason to make this stage reproducible afaikt
#slim -d "dir_out='abc_out/'" -d n_size=50 -d loci=200000 -d mu=3e-6 -d Na=600 -d N0=100 -d Nb=60 -d B_t=2 -d sfs1_shape=2 -d sfs1_mean=-2e-3 -d sfs2_shape=2 -d sfs2_mean=-2e-2 -d p_neutral=0.7 -d p_sfs1=0.2 src/nam_recap.slm

seeds = np.random.randint(0, int(2**62) - 1, n_draws)
#low informative prior set
Na = 1000 #ancestral pop size
#N0 = np.random.randint(Na, 10*Na, n_draws) #modern pop size
N0 = np.exp(np.random.uniform(np.log(Na), np.log(20*Na), n_draws)).astype(int)
Nb = np.random.randint(0.05*Na, Na, n_draws) #instant bottleneck pop size
B_t = int(0.067*Na)  #time after bottleneck
mut_props = np.random.dirichlet((3, 2, 1), n_draws) #proportion of mutation types
p_neutral = mut_props[:,0] #only need first two
p_sfs1 = mut_props[:,1]
sfs1_mean = -np.random.uniform(0, 0.01, n_draws).astype(np.half) #DFE negative only!
sfs1_shape = np.random.uniform(2, 100, n_draws).astype(np.half) #DFE positive only!
sfs2_mean = -np.random.uniform(0, 0.01, n_draws).astype(np.half) #DFE negative only!
sfs2_shape = np.random.uniform(2, 100, n_draws).astype(np.half) #DFE positive only!

#print(sfs2_shape)
stats_out = []

#params_file = open("abc_out/nam_abc_params.txt", "w")
#print(f"sample loci Na N0 Nb B_t B_r del_mu gamma_shape gamma_mean", file = params_file)
#seed=1616882522576_sample=50_loci=200000_mu=3e-06_Na=600_N0=100_Nb=60_Bt=2_sfs1shape=2_sfs1mean=-0.002_sfs2_shape=2_sfs2mean=-0.02_pneutral=0.7_psfs1=0.2_sumstats.txt

for i in range(n_draws):
    param_str = f"abc_out/seed__{seeds[i]}_Na__{Na}_N0__{N0[i]}_Nb__{Nb[i]}_Bt__{B_t}_sfs1shape__{sfs1_shape[i]}_sfs1mean__{sfs1_mean[i]}_sfs2shape__{sfs2_shape[i]}_sfs2mean__{sfs2_mean[i]}_pneutral__{p_neutral[i]}_psfs1__{p_sfs1[i]}_sumstats.txt"
    stats_out.append(param_str)

#slim -d "dir_out='abc_out/'" -d n_size=50 -d loci=200000 -d mu=3e-6 -d Na=600 -d N0=100 -d Nb=60 -d B_t=2 -d sfs1_shape=2 -d sfs1_mean=-2e-3 -d sfs2_shape=2 -d sfs2_mean=-2e-2 -d p_neutral=0.7 -d p_sfs1=0.2 src/nam_recap.slm
rule all:
    input:
        stats_out

rule run_slim:
    input:
        "src/nam_rawsfs.slim",
    output:
        "abc_out/seed__{seeds}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_pneutral__{p_neutral}_psfs1__{p_sfs1}_sumstats.txt"
    params:
        seeds = "{seeds}" ,
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{B_t}",
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}", 
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}", 
        pneutral="{p_neutral}", psfs1="{p_sfs1}"
    shell:
        """
        slim -s {params.seeds} -define rr={rr} \
        -define n_size={sample} -define loci={loci} -define mu={mu} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} \
        -define p_neutral={params.pneutral} -define p_sfs1={params.psfs1} {input} > {output}
        """
