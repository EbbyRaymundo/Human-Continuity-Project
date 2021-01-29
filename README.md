# Human-Continuity-Project

## Pipeline

    The following is the order in which you will typically use this pipeline.
1. Run PreProcessReads.py
    - Bam files need to be sorted and indexed
    - This file will take a very long time to run. Consider running this in the **background**
    - Use createAncientReads() to create fresh file of ancient individuals
    - AncientReads.output serves as a master file. createAncientReads() only ever needs to be used **once**
    - Use appendtoAncientReads() if AncientReads.output file exists and add a new column of reads for new ancient individuals
2. Run computeAlleleFreq.py
    - Use computeAlleleFrequency() with modern populations you want to select out of your .ind file
    - Use appendAncientIndividuals() to add a column for each ancient individual's reads from the AncientReads.output file
    - **If** the names of your bam files do not match those of your ind file, you can use a python dictionary in appendAncientIndividuals() with the bam file name (omitting the file extension) **first** and their name in the ind file **second**.
3. Run Schraiber's software. Here is a typical use example but the [Schraiber Github documentation](https://github.com/Schraiber/continuity/blob/master/README.md).
    - Output of this file is made using print statements (sorry, this was my best method for exporting the results). Edit this how you'd like.
```
# reading in data
unique_pops, inds, label, pops, freqs, read_lists = a_g.parse_reads_by_pop("reads/" + group + '_' + individual +".reads", "ind/" + group+ '_' + individual + ".ind")
# estimating parameters
opts_cont_false = a_g.optimize_pop_params_error_parallel(freqs,read_lists,48,continuity=False) #will only use one core; you can change the 1 to however many cores you want to use
print(opts_cont_false)
opts_cont_true = a_g.optimize_pop_params_error_parallel(freqs,read_lists,48,continuity=True)
print(opts_cont_true)
likelihood_false = np.array([-x[1] for x in opts_cont_false]) #minus sign is because scipy.optimize minimizes the negative log likelihood
likelihood_true = np.array([-x[1] for x in opts_cont_true])
LRT = 2*(likelihood_false - likelihood_true)
p_vals = scipy.stats.chi2.logsf(LRT,1) #returns the LOG p-values
```


```
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
```
4. Run CleanResults.py
    - This file has 2 modes. One creates a csv file with createCSV() and the other removes the "Reading line" text from your results but keeps everything else as a raw output with cleanResults(). Making any edits to the prints in the example above will require changing the createCSV() function.