import ancient_genotypes as a_g
import numpy as np
import scipy

unique_pops, inds, label, pops, freqs, read_lists = a_g.parse_reads_by_pop("ASW_HRR051935_HRR051936.reads", "testSchraiber.ind")

opts_cont_false = a_g.optimize_pop_params_error_parallel(freqs,read_lists,1,continuity=False)
opts_cont_true = a_g.optimize_pop_params_error_parallel(freqs,read_lists,1,continuity=True)
likelihood_false = np.array([-x[1] for x in opts_cont_false]) #minus sign is because scipy.optimize minimizes the negative log likelihood
likelihood_true = np.array([-x[1] for x in opts_cont_true])
LRT = 2*(likelihood_false - likelihood_true)
p_vals = scipy.stats.chi2.logsf(LRT,1) #returns the LOG p-values

print("Likelihood ratio test:", LRT)
print("P values (log(p values)):", p_vals)
