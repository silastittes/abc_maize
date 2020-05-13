import numpy as np
import random


#assume anc. Ne is 1000. median theta ~ 0.008. 4 * 1000 * 2e-6 = 0.008
#best estimate of per base mu is 3e-8. Going to use 2e-6
#best median estimate of recombination is r 1.6e-08 
#n = 3e-6/3e-8 = 67
#suggested rescaling is (1/2)*(1 - (1-2*r)^67) ~ 1.05e-6

#fixed params
loci = 20000000 #!!!
sample = 26
mu = 2e-6
rr = 1.05e-6

#priors
#n_draws = 150000
n_draws = 90000

#random.seed(214125) #random seed per sim is reported, so no reason to make this stage reproducible afaikt

seeds = np.random.randint(0, int(2**62) - 1, n_draws)
#low informative prior set
Na = 1000 #ancestral pop size
#N0 = np.random.randint(Na, 10*Na, n_draws) #modern pop size
N0 = np.exp(np.random.uniform(np.log(Na), np.log(20*Na), n_draws)).astype(int)
Nb = np.random.randint(0.05*Na, Na, n_draws) #instant bottleneck pop size
B_t = int(0.067*Na)  #time after bottleneck
mu_sv = np.random.uniform(0, 5e-8, n_draws)#.astype(np.half)
sfs1_mean = -np.random.uniform(0, 0.1, n_draws)#.astype(np.half) #DFE negative only!
sfs1_shape = np.random.uniform(2, 100, n_draws)#.astype(np.half) #DFE positive only!
sfs2_mean = -np.random.uniform(0, 0.1, n_draws)#.astype(np.half) #DFE negative only!
sfs2_shape = np.random.uniform(2, 100, n_draws)#.astype(np.half) #DFE positive only!


stats_out = []

for i in range(n_draws):
    param_str = f"abc_out/seed__{seeds[i]}_Na__{Na}_N0__{N0[i]}_Nb__{Nb[i]}_Bt__{B_t}_musv__{mu_sv[i]}_sfs1shape__{sfs1_shape[i]}_sfs1mean__{sfs1_mean[i]}_sfs2shape__{sfs2_shape[i]}_sfs2mean__{sfs2_mean[i]}_sumstats.txt"
    stats_out.append(param_str)

rule all:
    input:
        stats_out

rule run_slim:
    input:
        "src/nam_exons_rawsfs.slim",
    output:
        "abc_out/seed__{seeds}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt"
    params:
        seeds = "{seeds}" ,
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{B_t}",
        mu_sv = "{mu_sv}",
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}", 
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}" 
        #pneutral="{p_neutral}", psfs1="{p_sfs1}"
    shell:
        """
        slim -s {params.seeds} -define rr={rr} \
        -define n_size={sample} -define loci={loci} \
        -define mu={mu} -define mu_sv={params.mu_sv} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} \
        -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} {input} > {output}
        """
