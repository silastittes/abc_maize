import numpy as np
import random



#fixed params
loci = 2000000
sample = 50
mu = 3e-6
rr = 1.6e-7
#priors
n_draws = 10
random.seed(7141)

#slim -d "dir_out='abc_out/'" -d n_size=50 -d loci=200000 -d mu=3e-6 -d Na=600 -d N0=100 -d Nb=60 -d B_t=2 -d sfs1_shape=2 -d sfs1_mean=-2e-3 -d sfs2_shape=2 -d sfs2_mean=-2e-2 -d p_neutral=0.7 -d p_sfs1=0.2 src/nam_recap.slm

seeds = np.random.randint(0, int(2**62) - 1, n_draws)
#low informative prior set
Na = np.random.randint(100, 2000, n_draws) #ancestral pop size
N0 = np.random.randint(100, 2000, n_draws) #modern pop size
Nb = np.random.randint(10, 100, n_draws) #instant bottleneck pop size
B_t = np.random.randint(10, 100, n_draws)  #time after bottleneck, probably a bad prior
#mu = np.random.uniform(1e-8, 1e-6, n_draws).astype(np.half) #deleterious mutation rate
mut_props = np.random.dirichlet((2, 1, 1), n_draws) #proprotion of mutation types
p_neutral = mut_props[:,0] #only need first two
p_sfs1 = mut_props[:,1]
sfs1_mean = -np.random.uniform(0, 0.1, n_draws).astype(np.half) #DFE negative only!
sfs1_shape = np.random.randint(2, 10, n_draws).astype(np.half) #DFE positive only!
sfs2_mean = -np.random.uniform(0, 0.1, n_draws).astype(np.half) #DFE negative only!
sfs2_shape = np.random.randint(2, 10, n_draws).astype(np.half) #DFE positive only!

print(sfs2_shape)
stats_out = []

#params_file = open("abc_out/nam_abc_params.txt", "w")
#print(f"sample loci Na N0 Nb B_t B_r del_mu gamma_shape gamma_mean", file = params_file)
#seed=1616882522576_sample=50_loci=200000_mu=3e-06_Na=600_N0=100_Nb=60_Bt=2_sfs1shape=2_sfs1mean=-0.002_sfs2_shape=2_sfs2mean=-0.02_pneutral=0.7_psfs1=0.2_sumstats.txt

for i in range(n_draws):
    param_str = f"abc_out/seed={seeds[i]}_Na={Na[i]}_N0={N0[i]}_Nb={Nb[i]}_Bt={B_t[i]}_sfs1shape={sfs1_shape[i]}_sfs1mean={sfs1_mean[i]}_sfs2shape={sfs2_shape[i]}_sfs2mean={sfs2_mean[i]}_pneutral={p_neutral[i]}_psfs1={p_sfs1[i]}_sumstats.txt"
    stats_out.append(param_str)


#slim -d "dir_out='abc_out/'" -d n_size=50 -d loci=200000 -d mu=3e-6 -d Na=600 -d N0=100 -d Nb=60 -d B_t=2 -d sfs1_shape=2 -d sfs1_mean=-2e-3 -d sfs2_shape=2 -d sfs2_mean=-2e-2 -d p_neutral=0.7 -d p_sfs1=0.2 src/nam_recap.slm
rule all:
    input:
        stats_out

rule run_slim:
    input:
        "src/nam_recap.slm",
    output:
        "abc_out/seed={seeds}_Na={Na}_N0={N0}_Nb={Nb}_Bt={B_t}_sfs1shape={sfs1_shape}_sfs1mean={sfs1_mean}_sfs2shape={sfs2_shape}_sfs2mean={sfs2_mean}_pneutral={p_neutral}_psfs1={p_sfs1}_sumstats.txt"
    params:
        dir_out = "abc_out/",
        seeds = "{seeds}" ,
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{B_t}",
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}", 
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}", 
        pneutral="{p_neutral}", psfs1="{p_sfs1}"
    shell:
        """
        slim -s {params.seeds} -define 'dir_out="{params.dir_out}"' -define rr={rr} \
        -define n_size={sample} -define loci={loci} -define mu={mu} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} \
        -define p_neutral={params.pneutral} -define p_sfs1={params.psfs1} src/nam_recap.slm > {output}
        """
