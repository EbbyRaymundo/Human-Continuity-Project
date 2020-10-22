######################################################
#
# Author: Ebby Raymundo
# Date: 6/17/2020
# Description:
#   All files used are in eigenstrat format.
#   Individuals are searched v42 EG ind
#   file and their line number is used in
#   the v42 EG geno file as an index on each
#   line. Derived allele frequencies are
#   outputted as a tab delimited file 
#   (representing a table) in the line format:
#   "rsID | Chrom | Pos | Ancestral Allele | Frequency Derived Allele | Total individuals"
#   
#
######################################################

import numpy as np
import subprocess as sp

'''
Reads fileName and adds line number (starting at 0) of specified population to a list

fileName: ind file in eigenstrat format
group: search term. Uses .find() method. Is case sensitive.

return: list of line numbers for searched individuals
'''
def readIndividuals(fileName, group):

    file = open(fileName, 'r')
    lineNumber = 0
    individuals = [] # list that will hold line #'s of search individuals (will act as indexor later)

    for line in file:

        if (line.find(group) != -1): # search term is found in line
            individuals.append(lineNumber)  # add line # to list (will act as an index later)

        # advance
        lineNumber += 1

    file.close() # for good practice

    return(individuals)

'''
Uses string representing terminal command to run subprocess and return
the output as a list of tokens

command: String representing terminal command

return: list of string tokens from terminal output
'''
def runSubprocess(command):
    
    runCommand = sp.Popen(command.split(), stdout = sp.PIPE)
    terminalOutput = runCommand.communicate() # tuple of stdout (output, stderr)

    # Formats output into usable encoding and trims extra chars so we'll just have our
    # terminal output now. 
    stdoutList = terminalOutput[0].decode("utf-8").split('\n') # should have ["terminal output", ""] stored

    tokens = stdoutList[0].split() # tokenizes terminal output into list of strings

    return(tokens)

'''
Use line numbers from ind file to search geno file for allele freqs for a single SNP.

genoFile: geno file in eigenstrat format
indFile: ind file in eigenstrat format
snpFile: snp file in eigenstrat format
group: Search term to select populations from ind file. To specify no group, enter ''

return: numpy array of freqs, # of lines in geno file (# SNP's), searched group name
'''
def computeAlleleFreq(genoFile, indFile, snpFile, group):

    individualIndices = readIndividuals(indFile, group) # creates list of indices

    # need to count number of lines in .geno file to construct
    # the numpy array with the correct size. Appending to the
    # array requires resizing and becomes EXTREMELY expensive.

    lineCount = runSubprocess(f"wc -l {genoFile}") # contains ["lineCount", "genoFile"]

    gFile = open(genoFile, 'r')
    sFile = open(snpFile, 'r')
    
    outFile = open(f"{group}.output", 'w') # will overwrite contents of existing file
    outFile.write(f"Chrom\tPos\tAF\n") # creates header, aFileHeader[10:] skips chrom and pos

    freqs = np.zeros(int(lineCount[0])) # one element for each SNP freq (line)
    lineNumber  = 0  # represents which SNP we're on
    snpInfoLine = sFile.readline() # for loop iterates over geno file, iterate over SNP file simultaneously
    
    # for each line in file, calculate freq of SNP
    for line in gFile:

        alleleTotal = 0

        # for each char in line, a 0 means no copies of reference allele, 1 is one
        # copy, 2 is two copies, and 9 is no data. 9's are handled by excluding from
        # total allele count.
        for index in individualIndices:

            if line[index] != '9': # if read has data

                freqs[lineNumber] += float(line[index]) # need to cast char to float for numpy array
                alleleTotal += 2 # we're diploid, two copies of each allele

        freqs[lineNumber] = 1 - (freqs[lineNumber] / alleleTotal) # calculate derived allele freq step

        # need to skip alleles at frequency 0 or 1 in reference modern
        # population for Schraiber's program. Skips writing out line
        # for SNP and moves to next iteration. Also skipping if snp
        # has no data for calculating AF
        if (freqs[lineNumber] == 0 or freqs[lineNumber] == 1 or np.isnan(freqs[lineNumber])):
            # need to advance these before passing onto next iteration
            snpInfoLine = sFile.readline()
            lineNumber += 1 # advance to next SNP's counts
            continue

        lineList = snpInfoLine.split() # should be in "rsID | chromosome | pos | physPos | refAllele | newAllele" format
        
        # writing chrom, physpos, and allele freq
        outFile.write(f"{lineList[1]}\t{lineList[3]}\t{freqs[lineNumber]}\n")
        snpInfoLine = sFile.readline() # advance to info of next SNP

        lineNumber += 1 # advance to next SNP's counts
            
    gFile.close()
    sFile.close() # for good practice
    outFile.close()

    return(freqs, lineNumber, group)


'''
Appends the read data for an ancient individual to the
1k genome population data whose allele frequencies
have already been calculated.

group: 1k genomes group to append to
individuals: List of names of ancient individuals. Name should 
             match those in ind file of ancient samples.
ancientFile: Preprocessed ancient individual data (from running 
             PreProcessReads.py). Contains mpileup read counts 
             for each snp in snpFile.
'''
def appendAncientIndividual(group, individuals, ancientFile):

    # contains header "Chrom | Pos | anc1_der | anc1_anc | anc1_other | anc2_der | ..."
    aFile = open(ancientFile, 'r')
    aFileLine = aFile.readline().split() # sitting on first line, is now in tokens
    indList = [] # for keeping track of where individuals we want are in aFileLine list

    for i in range(len(individuals)):
        indList.append(aFileLine.index(f"{individuals[i]}_der")) # adds first index of ancient individual

    groupFile = open(f"{group}.output", 'r') # opens 1k genome source file
    outFile = group

    for i in range(len(individuals)):
        outFile = f"{outFile}_{individuals[i]}" # will finish as "group_ind1_ind2_..." for file name

    outFile = open(f"{outFile}.reads", "w")

    header = groupFile.readline() # won't start on header line in main loop
    header = header.strip()

    for i in range(len(indList)):
        header = f"{header}\t{aFileLine[indList[i]]}\t{aFileLine[indList[i] + 1]}\t{aFileLine[indList[i] + 2]}" # adds anci_der anci_anc anci_other

    outFile.write(f"{header}\n")

    aFileLine = aFile.readline() # sitting on first data line
    

    # for each line of 1k genome group, advance ancient reads file until
    # on correct allele. We skipped alleles in cases where allele was
    # was fixed, extinct, or had no data.
    for line in groupFile:

        lineList = line.split() # in format "chrom | pos | AF"
        aFileLine = aFileLine.split() # in format "chrom | pos | anc1_der | anc1_anc | anc1_other | anc2_der..."

        # if chrom or allele position don't match, need to advance
        # to next line of ancient reads until they both do

        while (lineList[0] != aFileLine[0] or lineList[1] != aFileLine[1]):
            aFileLine = aFile.readline()
            aFileLine = aFileLine.split()

        writtenLine = f"{lineList[0]}\t{lineList[1]}\t{lineList[2]}"

        for i in range(len(indList)):
            writtenLine = f"{writtenLine}\t{aFileLine[indList[i]]}\t{aFileLine[indList[i] + 1]}\t{aFileLine[indList[i] + 2]}"

        outFile.write(f"{writtenLine}\n")
        aFileLine = aFile.readline()

    aFile.close()   # for good practice
    outFile.close()

    return(indList)


#####################################################
#
# Main
#
#####################################################

#eigenstratIndFile   = "v42.4.1240K.EG.ind"
#eigenstratGenoFile  = "v42.4.1240K.EG.geno"


'''
testIndFile = "test.ind"
testGenoFile = "test.geno"
testSNPFile = "test.snp"
searchTerm = "test"
chimpFile = "testChimp.geno"
reads = "AncientReads.output"
computeAlleleFreq(testGenoFile, testIndFile, testSNPFile, searchTerm, reads)
'''

IndFile = "v42.4.1240K.EG.ind"
GenoFile = "v42.4.1240K.EG.geno"
SNPFile = "v42.4.1240K.EG.snp"
searchTerm = "ASW"
reads = "AncientReads.output"

#searchTerms = ["ACB", "ASW","BEB", "GBR", "CDX", "CLM", "ESN", "FIN", "GWD", "GIH", "CHB", "CHS", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PUR", "PJL", "STU", "TSI", "YRI", "CEU"]
#computeAlleleFreq(GenoFile, IndFile, SNPFile, searchTerm)
print(appendAncientIndividual("ASW", ["HRR051935", "HRR051936"], "AncientReads.output"))
#for i in range(len(searchTerms)):
#    computeAlleleFreq(GenoFile, IndFile, SNPFile, searchTerms[i])