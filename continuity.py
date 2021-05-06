import ancient_genotypes as a_g
import scipy.stats
import numpy as np

modern_pops = ["CHB", "CHS", "CDX", "JPT", "KHV", "CEU"]
'''
ancient_Individuals = ["HRR051935", "HRR051936", "HRR051937",\
                      "HRR051938", "HRR051939", "HRR051940",\
                      "HRR051941", "HRR051942", "HRR051943",\
                      "HRR051944", "HRR051945", "HRR051946",\
                      "HRR051947", "HRR051948", "HRR051949",\
                      "HRR051950", "HRR051951", "HRR051952",\
                      "HRR051954", "HRR051955",\
                      "HRR051956", "HRR051958",\
                      "HRR051959", "HRR051960"]
'''
'''
ancient_Individuals = ["HRR051935", "HRR051936", "HRR051937",
                        "HRR051938_HRR051939_HRR051940", "HRR051941",
                        "HRR051942", "HRR051943_HRR051944",
                        "HRR051945_HRR051946",
                        "HRR051947_HRR051948_HRR051949_HRR051950_HRR051951_HRR051952_HRR051954",
                        "HRR051955_HRR051956_HRR051958_HRR051959",
                        "HRR051960"]
'''
ancient_Individuals = ["HRR051935", "HRR051936", "HRR051937",
                        "HRR051938_HRR051939_HRR051940", "HRR051941",
                        "HRR051942", "HRR051943_HRR051944",
                        "HRR051945",
                        "HRR051947_HRR051948_HRR051949_HRR051950",
                        "HRR051955_HRR051956_HRR051958",
                        "HRR051960"]
'''
modern_pops = ["KHV", "CEU"]
ancient_Individuals = ["HRR051935"]
'''
for group in modern_pops:
    for individual in ancient_Individuals:

        # reading in data
        unique_pops, inds, label, pops, freqs, read_lists = a_g.parse_reads_by_pop("reads/" + group + '_' + individual +".reads", "ind/" + group + '_' + individual + ".ind")

        a_g.coverage_filter(read_lists) # removing extremely high and low coverage SNPs

        # estimating parameters
        opts_cont_false = a_g.optimize_pop_params_error_parallel(freqs,read_lists,48,continuity=False) #will only use one core; you can change the 1 to however many cores you want to use
        print(opts_cont_false)
        opts_cont_true = a_g.optimize_pop_params_error_parallel(freqs,read_lists,48,continuity=True)
        print(opts_cont_true)
        likelihood_false = np.array([-x[1] for x in opts_cont_false]) #minus sign is because scipy.optimize minimizes the negative log likelihood
        likelihood_true = np.array([-x[1] for x in opts_cont_true])
        LRT = 2*(likelihood_false - likelihood_true)
        p_vals = scipy.stats.chi2.logsf(LRT,1) #returns the LOG p-values

        print("1k genomes group: " + group)
        print("Ancient Individual: " + individual)
        print("t1 continuity false:")
        print(opts_cont_false[0][0])
        print("t2 continuity false: ")
        print(opts_cont_false[0][1])
        print("continuity false error: ")
        for error in opts_cont_false[2:]:
            print(error)
        print("t1 continuity true: ")
        print( opts_cont_true[0][0])
        print("t2 continuity true:")
        print("0")
        print("continuity true error: ")
        for error in opts_cont_true[1:]:
            print(error)
        print("LRT: ")
        print(LRT)
        print("P values: ")
        print(p_vals)
        print('')

