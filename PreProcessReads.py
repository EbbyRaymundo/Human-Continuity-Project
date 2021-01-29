########################################
#
# Author: Ebby Raymundo
# Date: 9/23/20
# Description:
#   This program preprocesses the reads 
#   on all ancient individuals provided.
#   This is so that the time consuming
#   mpileup step can be done all in one
#   run. The file can be read later and
#   have data pulled from necessary
#   alleles. Requires corresponding
#   eigenstrat SNP file
#
########################################
import subprocess as sp
import pysam

'''
Uses string representing mpileup command to run subprocess and return
the output as a list of tokens

command: String representing mpileup command

return: string list of tokens from terminal output
'''
def runSubprocess(command):

    runCommand = sp.Popen(command.split(), stdout = sp.PIPE)
    terminalOutput = runCommand.communicate() # tuple of stdout (output, stderr)

    # Formats output into usable encoding and trims extra chars so we'll just have our
    # terminal output now. 
    stdoutList = terminalOutput[0].decode("utf-8").split('\n') # should have ["terminal output", ""] stored

    tokens = stdoutList[0].split('\t') # now a list of strings ["chrom", "pos", "N", "total reads", "reads", "read quality"]

    return (tokens)

'''
Creates new file containing reads of each ancient
sample's mpileup reads. WILL OVERWRITE EXISTING
FILE IF RUN BEFORE. In format:

chrom | pos | anc1_der | anc1_anc | anc1_other | anc2_der | ...

individuals: list of ancient individuals who have sorted and
             indexed bam files available
snpFile: eigenstrat format snp file
bamFilePath: directory path to the bam files. Make sure that
             the path you give it is the FOLDER, not a file.
'''

def createAncientReads(individuals, snpFile, bamFilePath):

    sFile = open(snpFile, 'r')
    outFile = open("AncientReads.output", 'w') # will OVERWRITE contents of existing file

    header = 'Chrom\tPos'

    for i in range(len(individuals)): # need to create header for all individuals in list

        header = f"{header}\t{individuals[i]}_der\t{individuals[i]}_anc\t{individuals[i]}_other"

    outFile.write(f"{header}\n")

    for line in sFile:

        lineList = line.split() # should be in "rsID | chrom | pos | physPos | refAllele | newAllele" format
        writtenLine = f"{lineList[1]}\t{lineList[3]}" # resets/starts the line we'll write to AncientReads.output file

        for i in range(len(individuals)): # run mpileup on each individual

            # in format "samtools mpileup -r chromosome:pos-pos bamFile"
            # Even if you only want one positon, do pos-pos. Otherwise
            # You'll get reads from that position to the end of the
            # chromosome.

            if (lineList[1] == '23'): # if X chromosome
                reads = pysam.mpileup('-r', f"X:{lineList[3]}-{lineList[3]}", f"{bamFilePath}{individuals[i]}.sorted.bam")
                reads = reads.split()

            elif(lineList[1] == '24'): # if Y chromosome
                reads = pysam.mpileup('-r', f"Y:{lineList[3]}-{lineList[3]}", f"{bamFilePath}{individuals[i]}.sorted.bam")
                reads = reads.split()

            else:
                reads = pysam.mpileup('-r', f"{lineList[1]}:{lineList[3]}-{lineList[3]}", f"{bamFilePath}{individuals[i]}.sorted.bam")
                reads = reads.split()
            # list of strings ["chrom", "pos", "N", "total reads", "reads", "read quality"]

            if not reads: # no read was found
                ancReads = 0
                derReads = 0
                otherReads = 0

            else:
                totalReads = int(reads[3])
                ancReads = reads[4].count(lineList[4])
                derReads = reads[4].count(lineList[5])
                otherReads = totalReads - ancReads - derReads

            writtenLine = f"{writtenLine}\t{derReads}\t{ancReads}\t{otherReads}" # append reads to line

        outFile.write(f"{writtenLine}\n")

    sFile.close()   # for good practice
    outFile.close()

########################################
#
# Main
#
########################################

# test bois
#ancientIndividuals = ["HRR051935"]
#snpFile = "test.snp"

snpFile = "v42.4.1240K.EG.snp"

ancientIndividuals = ["HRR051935", "HRR051936", "HRR051937",\
                      "HRR051938", "HRR051939", "HRR051940",\
                      "HRR051941", "HRR051942", "HRR051943",\
                      "HRR051944", "HRR051945", "HRR051946",\
                      "HRR051947", "HRR051948", "HRR051949",\
                      "HRR051950", "HRR051951", "HRR051952",\
                      "HRR051953", "HRR051954", "HRR051955",\
                      "HRR051956", "HRR051957", "HRR051958",\
                      "HRR051959", "HRR051960"]

bamFiles = "/home/classes/myanglab/data/earlyCN/"

createAncientReads(ancientIndividuals, snpFile, bamFiles)
