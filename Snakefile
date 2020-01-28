import numpy as np
import random
from scipy.spatial import distance
import msprime, pyslim
import tskit

#functions

#need to pass more info from file names to function!!!
#Ne, N0, Na,
#sn = "{sn}", l = "{l}", a = "{a}", o = "{o}",b = "{b}", t = "{t}",r = "{r}", u= "{u}", s= "{s}", gm= "{gm}"
#sample{sn}_loci{l}_Na{a}_N0{o}_Nb{b}_Bt{t}_Br{r}_del_mu{u}_gammaShape{s}_gammaMean{gm}_sumstat.txt 
def recap_sumstats(tree_file, sum_stat_file, 
    recombination_rate, sample,
    loci, Na, N0, Nb, Bt, Br, mu, del_mu, gammaShape, gammaMean):

    sfile = open(sum_stat_file, 'w')

    slim_ts = pyslim.load(tree_file)

    recap = slim_ts.recapitate(recombination_rate=recombination_rate, Ne=Na)

    ts = pyslim.SlimTreeSequence(
          msprime.mutate(recap, rate=mu, keep=True))

    sample_pop = ts.samples()[(np.random.choice(range(N0), sample, replace = False))]

    sfs = ts.allele_frequency_spectrum([sample_pop], polarised=False,  span_normalise=False)
    sfs = sfs/sum(sfs)
    sfs = sfs[0:int((len(sfs) + 1*len(sfs)%2)/2)]

    sfs_cat = list(range(len(sfs)))
    
    sfs_pairs = distance.pdist(np.reshape(np.multiply(sfs, sfs_cat), (len(sfs),1)), metric = "euclidean")

    td = tskit.TreeSequence.Tajimas_D(ts, sample_sets = [sample_pop])[0]
    pi = ts.pairwise_diversity(sample_pop) / loci #0.0053161183673468655
    sfs_dist_mean = sfs_pairs.mean()
    sfs_dist_sd = sfs_pairs.std()
   
    print(f"{Na} {N0} {Nb} {Bt} {Br} {del_mu} {gammaShape} {gammaMean} {td} {pi} {sfs_dist_mean} {sfs_dist_sd}",
          file = sfile)

####RANDOM SEED. FIXES OUTPUT FILES

#to do
    #get confifuring to NOT redo analysis every time?

#fixed params
loci = 1000000
true_mu = 3e-8
true_r = 1.6e-08
mu = 3e-7
sc_fct = mu/true_mu
recombination_rate=(1/2)*(1 - (1-2*true_r)**sc_fct) 
sample = 50


#priors
n_draws = 10
random.seed(71415)
seeds = np.random.randint(0, int(2**62) - 1, n_draws)
#low informative prior set
Na = np.random.randint(10, 10000, n_draws) #ancestral pop size
N0 = np.random.randint(50, 10000, n_draws) #modern pop size
Nb = np.random.randint(2, 10, n_draws) #instant bottleneck pop size
B_t = np.random.randint(100, 1000, n_draws)  #time after bottleneck, probably a bad prior
B_r = np.random.randint(10, 1000, n_draws) #time after bottleneck, probably a bad prior
del_mu = np.random.uniform(1e-9, 1e-6, n_draws) #deleterious mutation rate
print(del_mu)
gamma_mean = -np.random.uniform(0, 0.1, n_draws) #DFE negative only!
gamma_shape = np.random.randint(2, 20, n_draws) #DFE positive only!

sims_out = []
stats_out = []

#params_file = open("abc_out/nam_abc_params.txt", "w")
#print(f"sample loci Na N0 Nb B_t B_r del_mu gamma_shape gamma_mean", file = params_file)

for i in range(n_draws):
    param_str = f"{seeds[i]} {sample} {loci} {Na[i]} {N0[i]} {Nb[i]} {B_t[i]} {B_r[i]} {del_mu[i]} {gamma_shape[i]} {gamma_mean[i]}"
    #print(param_str, file = params_file)

    c_file = f"abc_out/seed{seeds[i]}_sample{sample}_loci{loci}_Na{Na[i]}_N0{N0[i]}_Nb{Nb[i]}_Bt{B_t[i]}_Br{B_r[i]}_del_mu{del_mu[i]}_gammaShape{gamma_shape[i]}_gammaMean{gamma_mean[i]}"
    sims_out.append(f"{c_file}.trees")
    stats_out.append(f"{c_file}_sumstat.txt")

rule all:
    input:
        sims_out,
        stats_out

rule run_slim:
    input:
        "src/nam_recap.slm",
    output:
        "abc_out/seed{sd}_sample{sn}_loci{l}_Na{a}_N0{o}_Nb{b}_Bt{t}_Br{r}_del_mu{u}_gammaShape{s}_gammaMean{gm}.trees"
    params:
        sd = "{sd}", sn = "{sn}", l = "{l}", a = "{a}", o = "{o}",b = "{b}", t = "{t}",r = "{r}", u= "{u}", s= "{s}", gm= "{gm}"
    shell:
        """
slim -seed {params.sd} -d sample={params.sn} -d loci={params.l} -d Na={params.a} -d N0={params.o} -d Nb={params.b} -d B_t={params.t} -d B_r={params.r} -d mu={params.u} -d gamma_shape={params.s} -d gamma_mean={params.gm} -d 'tree_file="{output}"' {input}
        """


rule sumstat:
    input:
        infile = "abc_out/seed{sd}_sample{sn}_loci{l}_Na{a}_N0{o}_Nb{b}_Bt{t}_Br{r}_del_mu{u}_gammaShape{s}_gammaMean{gm}.trees"
    output:
        outfile = "abc_out/seed{sd}_sample{sn}_loci{l}_Na{a}_N0{o}_Nb{b}_Bt{t}_Br{r}_del_mu{u}_gammaShape{s}_gammaMean{gm}_sumstat.txt"
    params:
        sn = "{sn}", l = "{l}", a = "{a}", o = "{o}",b = "{b}", t = "{t}",r = "{r}", u= "{u}", s= "{s}", gm= "{gm}"
    run:
        recap_sumstats(input.infile, output.outfile, 
                        recombination_rate = recombination_rate, sample = sample, loci = loci, 
                        Na = int(params.a), N0 = int(params.o), Nb = int(params.b),
                        Bt = int(params.t), Br = int(params.r), mu = mu, 
                        del_mu = float(params.u), gammaShape = float(params.s), gammaMean = float(params.gm))

