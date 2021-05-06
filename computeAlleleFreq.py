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

return: list of line numbers for searched individuals,
        list of lines that matched search

'''
def readIndividuals(fileName, group):

    file = open(fileName, 'r')
    lineNumber = 0
    individuals = [] # list that will hold line #'s of search individuals (will act as indexor later)
    lineList = [] # to keep track of what we're searching. Can serve as debug tool, but is actually
                  # used in appendAncientIndividuals() function

    for line in file:

        if group in line:
            individuals.append(lineNumber)  # add line # to list (will act as an index later)
            lineList.append(line)

        # advance
        lineNumber += 1

    file.close() # for good practice

    return(individuals, lineList)

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

    individualIndices, searchedLines = readIndividuals(indFile, group) # creates list of indices

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
    snpLost = 0 # keeps track of how many SNP's had no data for the group
    
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

            else:
                snpLost += 1

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

    return(freqs, lineNumber, group, snpLost)


'''
Appends the read data for an ancient individual to the
1k genome population data whose allele frequencies
have already been calculated using computeAlleleFrequency().
File name will be a compound of 1k genomes group and
individuals using their bam file name. Also creates
a corresponding .ind file for the newly made reads
file (uses same naming convention).
Ex: group_individual1_individual2.reads
    group_individual1_individual2.ind

INPUTS
group: 1k genomes group to append to
individuals: LIST of names of ancient individuals. Names should 
             match those in header of ancientFile.
ancientFile: Preprocessed ancient individual data (from running 
             PreProcessReads.py). Contains mpileup read counts 
             for each snp in snpFile.
nameDict: Optional argument if your .bam file names don't match
          the .ind file for your ancient individuals. This will
          make the .reads header match your .ind file.
          Note: The .reads and .ind file name isn't impacted.
indFile: If you give a nameDict, you need to provide the
                corresponding .ind file.

OUTPUTS
.reads (outFile): has header group_ind1_ind2_... and contains
                  derived allele frequencies and read counts
                  of the snp from bam files
.ind (newIndFile): .ind file of desired alternative names for
                   individuals of the .bam files. Underscores
                   are removed from names within the file.
'''
def appendAncientIndividual(group, individuals, ancientFile, nameDict = None, indFile = None):

    # if we're given a nameDict, it's implied that the bamfile names
    # don't match what's in the .ind file for the ancient individual
    if (nameDict != None):
        if (indFile == None):
            print("You must provide a corresponding .ind file")
            return()

        newIndFile = group

        for each in individuals:
            newIndFile = f"{newIndFile}_{each}" # will finish as "group_ind1_ind2_..." for file name

        newIndFile = open(f"{newIndFile}.ind", "w")

        for each in individuals:
            print(each)
            print(nameDict[each])
            print(individuals)
            individualIndices, searchedLines = readIndividuals(indFile, nameDict[each])
            fixedLine = searchedLines[0].replace('_', '') # need to remove underscores
            newIndFile.write(fixedLine) # searchedLines should only contain 1 value

        newIndFile.close()

    # contains header "Chrom | Pos | anc1_der | anc1_anc | anc1_other | anc2_der | ..."
    aFile = open(ancientFile, 'r')
    aFileLine = aFile.readline().split() # sitting on first line, is now in tokens
    indList = [] # for keeping track of where individuals we want are in aFileLine list

    for each in individuals:
        indList.append(aFileLine.index(f"{each}_der")) # adds first index of ancient individual

    groupFile = open(f"{group}.output", 'r') # opens 1k genome source file
    outFile = group

    for each in individuals:
        outFile = f"{outFile}_{each}" # will finish as "group_ind1_ind2_..." for file name

    outFile = open(f"{outFile}.reads", "w")

    header = groupFile.readline() # won't start on header line in main loop
    header = header.strip()

    if (nameDict != None): # need to base header off dictionary and remove '_'

        fixedDict = {} # making a new dictionary keeps the argument intact

        for key in nameDict:
            fixedDict[key] = nameDict[key].replace('_', '')

        for i in range(len(indList)):
            # adds anci_der anci_anc anci_other
            # aFileLine[indList[i]] contains "individual_der".
            # Here we replace "individual" with name in the
            # dictionary but keep the "_der"
            header = f"{header}\t{aFileLine[indList[i]].replace(individuals[i], fixedDict[individuals[i]])}\t{aFileLine[indList[i] + 1].replace(individuals[i], fixedDict[individuals[i]])}\t{aFileLine[indList[i] + 2].replace(individuals[i], fixedDict[individuals[i]])}"

    else:

        for i in range(len(indList)):
            # adds anci_der anci_anc anci_other
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
#searchTerm = "CHB"
reads = "AncientReads.output"

#searchTerms = ["ACB", "ASW","BEB", "GBR", "CDX", "CLM", "ESN", "FIN", "GWD", "GIH", "CHB", "CHS", "IBS", "ITU", "JPT", "KHV", "LWK", "MSL", "MXL", "PEL", "PUR", "PJL", "STU", "TSI", "YRI", "CEU"]
searchTerms = ["CHB", "CHS", "CDX", "JPT", "KHV", "CEU"]

'''
for group in searchTerms:
     results = computeAlleleFreq(GenoFile, IndFile, SNPFile, group)
     print(group, results[3])
'''

ancientIndividuals = ["HRR051935", "HRR051936", "HRR051937",\
                      "HRR051938", "HRR051939", "HRR051940",\
                      "HRR051941", "HRR051942", "HRR051943",\
                      "HRR051944", "HRR051945", "HRR051946",\
                      "HRR051947", "HRR051948", "HRR051949",\
                      "HRR051950", "HRR051951", "HRR051952",\
                      "HRR051954", "HRR051955",\
                      "HRR051956", "HRR051958",\
                      "HRR051959", "HRR051960"]

ancientDict = {"HRR051935":"Yumin",\
               "HRR051936":"Bianbian",\
               "HRR051937":"BS",\
               "HRR051938":"XJS1309_M7",\
               "HRR051939":"XJS1311_M16",\
               "HRR051940":"XJS1309_M4",\
               "HRR051941":"Xiaogao",\
               "HRR051942":"Qihe2_d",\
               "HRR051943":"LD1",\
               "HRR051944":"LD2",\
               "HRR051945":"SuogangB1_d",\
               "HRR051946":"SuogangB3_d",\
               "HRR051947":"L5705",\
               "HRR051948":"L5700",\
               "HRR051949":"L5692_d",\
               "HRR051950":"L5706_d",\
               "HRR051951":"L5704_d",\
               "HRR051952":"L5703_d",\
               "HRR051954":"L5701_d",\
               "HRR051955":"L7415",\
               "HRR051956":"L7417_d",\
               "HRR051958":"L5698_d",\
               "HRR051959":"L5696_d",\
               "HRR051960":"L5694"}
'''
# removed HRR051957 and HRR051953 for not having a corresponding name in the ind file
for group in searchTerms:
    for each in ancientIndividuals:
        appendAncientIndividual(group, [each], reads, ancientDict, "ind/early_CN.ind")
'''

for group in searchTerms:
    appendAncientIndividual(group, ["HRR051935"], "AncientReads.output", ancientDict, "early_CN.ind") # yumin
    appendAncientIndividual(group, ["HRR051936"], "AncientReads.output", ancientDict, "early_CN.ind") # bianbian
    appendAncientIndividual(group, ["HRR051937"], "AncientReads.output", ancientDict, "early_CN.ind") # boshan
    appendAncientIndividual(group, ["HRR051938", "HRR051939", "HRR051940"], "AncientReads.output", ancientDict, "early_CN.ind") # Xiaojinshan
    appendAncientIndividual(group, ["HRR051941"], "AncientReads.output", ancientDict, "early_CN.ind") # XIaogao
    appendAncientIndividual(group, ["HRR051942"], "AncientReads.output", ancientDict, "early_CN.ind") # Qihe
    appendAncientIndividual(group, ["HRR051943", "HRR051944"], "AncientReads.output", ancientDict, "early_CN.ind") # Liangdao
    appendAncientIndividual(group, ["HRR051945"], "AncientReads.output", ancientDict, "early_CN.ind") # Suogang
    appendAncientIndividual(group, ["HRR051947", "HRR051948", "HRR051949", "HRR051950"], "AncientReads.output", ancientDict, "early_CN.ind") # Xitoucun
    appendAncientIndividual(group, ["HRR051955", "HRR051956", "HRR051958"], "AncientReads.output", ancientDict, "early_CN.ind") # Tanshishan
    appendAncientIndividual(group, ["HRR051960"], "AncientReads.output", ancientDict, "early_CN.ind") # Chuanyun
