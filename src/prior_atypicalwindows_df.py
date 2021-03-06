import numpy as np
import pandas as pd
import random


"""
BUILD COMMAND script
"""


np.random.seed(3683)

n_draws = 300
#n_draws = 10

#assume anc. Ne is 1000. median theta ~ 0.008. 4 * 1000 * 2e-6 = 0.008
#best estimate of per base mu is 3e-8. Going to use 2e-6
#best median estimate of recombination is r 1.6e-08 
#n = 3e-6/3e-8 = 67
#suggested rescaling is (1/2)*(1 - (1-2*r)^67) ~ 1.05e-6

Na = 1000
ns0 = 100
ngenes_small = np.random.randint(3, 300, int(n_draws/2))
ngenes_big = np.random.randint(300, 3000, int(n_draws/2))
ngenes_median = [300] * ns0
ngenes = np.concatenate((ngenes_small, ngenes_big, ngenes_median))

n_draws = len(ngenes)

print(ngenes)

prior_df = pd.DataFrame({
    #fixed params
    'ngenes' : ngenes,
    'loci' : [20000000] * n_draws,
    'sample' : [26] * n_draws,
    'mu' : [2e-6] * n_draws,
    'rr' : [1.05e-6] * n_draws,
    'seeds' : np.random.randint(0, int(2**62) - 1, n_draws),
    'Na' : [Na] * n_draws, #ancestral pop size
    'N0' : np.exp(np.random.uniform(np.log(Na), np.log(20*Na), n_draws)).astype(int),
    'Nb' : np.random.randint(0.01*Na, Na, n_draws), #instant bottleneck pop size
    'B_t' : int(0.067*Na),  #time after bottleneck
    'mu_sv' : np.random.uniform(0, 5e-8, n_draws), 
    'sfs1_mean' : np.concatenate((-np.random.uniform(0, 0.1, n_draws - ns0), [0] * ns0)),
    'sfs1_shape' : np.random.uniform(0, 100, n_draws),
    'sfs2_mean' :  np.concatenate((-np.random.uniform(0, 0.1, n_draws - ns0), [0] * ns0)),
    'sfs2_shape' : np.random.uniform(0, 100, n_draws)
})

prior_df.to_csv('prior_atypical_df.csv')
