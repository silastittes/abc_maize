import numpy as np
import pandas as pd
import random

#TO TRIGGER NEW RUN OF PIPELINE, RUN:
#$ python src/prior_atypicalwindows_df.py

prior_df = pd.read_csv('prior_atypical_df.csv').astype({'ngenes':'int64', 'Na': 'int64', 'N0': 'int64', 'Nb': 'int64'})

stats_out = [ f"atypical_abc_out/seed__{int(row['seeds'])}_ngenes__{int(row['ngenes'])}_Na__{int(row['Na'])}_N0__{int(row['N0'])}_Nb__{int(row['Nb'])}_Bt__{int(row['B_t'])}_musv__{row['mu_sv']}_sfs1shape__{row['sfs1_shape']}_sfs1mean__{row['sfs1_mean']}_sfs2shape__{row['sfs2_shape']}_sfs2mean__{row['sfs2_mean']}_sumstats.txt" for index, row in prior_df.iterrows()]

#fixed params
loci = list(prior_df['loci'])[0]
sample =  list(prior_df['sample'])[0]
mu =  list(prior_df['mu'])[0]
rr =  list(prior_df['rr'])[0]


#for i in range(n_draws):
#    param_str = f"abc_out/seed__{seeds[i]}_Na__{Na}_N0__{N0[i]}_Nb__{Nb[i]}_Bt__{B_t}_musv__{mu_sv[i]}_sfs1shape__{sfs1_shape[i]}_sfs1mean__{sfs1_mean[i]}_sfs2shape__{sfs2_shape[i]}_sfs2mean__{sfs2_mean[i]}_sumstats.txt"
#    stats_out.append(param_str)

rule all:
    input:
        stats_out #,"atypical_pop_sfs.txt", "atypical_pop.txt"

rule val_slim:
    input:
        slim = "src/nam_exonsSIZE_rawsfs.slim",
    output:
        slim_out = "atypical_abc_out/seed__{seeds}_ngenes__{ngenes}_Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt"
    params:
        seeds = "{seeds}" ,
        ngenes = "{ngenes}",
        Na = "{Na}", N0 = "{N0}",Nb = "{Nb}", Bt = "{B_t}",
        mu_sv = "{mu_sv}",
        sfs1shape="{sfs1_shape}", sfs1mean="{sfs1_mean}", 
        sfs2shape="{sfs2_shape}", sfs2mean="{sfs2_mean}" 
        #pneutral="{p_neutral}", psfs1="{p_sfs1}"
    shell:
        """
        slim -s {params.seeds} -define ngenes={params.ngenes} -define rr={rr} \
        -define n_size={sample} -define loci={loci} \
        -define mu={mu} -define mu_sv={params.mu_sv} \
        -define Na={params.Na} -define N0={params.N0} -define Nb={params.Nb} -define B_t={params.Bt} \
        -define sfs1_shape={params.sfs1shape} -define sfs1_mean={params.sfs1mean} \
        -define sfs2_shape={params.sfs2shape} -define sfs2_mean={params.sfs2mean} {input.slim} > {output.slim_out}
        """



#rule val_sim_df:
#    input:
#        expand("atypical_abc_out/seed__{seeds}_ngenes__{ngenes}__Na__{Na}_N0__{N0}_Nb__{Nb}_Bt__{B_t}_musv__{mu_sv}_sfs1shape__{sfs1_shape}_sfs1mean__{sfs1_mean}_sfs2shape__{sfs2_shape}_sfs2mean__{sfs2_mean}_sumstats.txt", 
#                        zip,
#                        seeds =  prior_df['seeds'],
#                        ngenes = prior_df['ngenes'],
#                        Na = prior_df['Na'],
#                        N0 = prior_df['N0'],
#                        Nb = prior_df['Nb'],
#                        B_t = prior_df['B_t'],
#                        mu_sv = prior_df['mu_sv'],
#                        sfs1_mean = prior_df['sfs1_mean'],
#                        sfs1_shape = prior_df['sfs1_shape'],
#                        sfs2_mean = prior_df['sfs2_mean'],
#                        sfs2_shape = prior_df['sfs2_shape'],
#                        )
#    output:
#        sfs = "atypical_pop_sfs.txt",
#        param = "atypical_pop.txt"
#    shell:
#        """
#            grep -h -B2 "SFS:" {input} | grep "SFS:" | cut -d ":" -f2- | sed 's/^ //g' > {output.sfs}
#            grep -h -B2 "SFS:" {input} | grep -h -A1 "#Na" | grep -v "\-\-" | grep -v "#Na" > {output.param} 
#        """

