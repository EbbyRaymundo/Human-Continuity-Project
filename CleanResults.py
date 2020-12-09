##################################
#
# Author: Ebby Raymundo
# Date: 12/4/2020
# Desc.: 
#   This file is used to remove
#   "Reading line:" segments from
#   your results file. This step
#   is unavoidable at the moment
#   because results can't be
#   written to a file without
#   encoding errors so output needs
#   to be piped into a file.
#
##################################
'''
Uses .results file to create a .csv that
contains the 1k genome population, ancient
individual ID, LRT (likelihood ratio test)
value, and the p value.
"population" | "individual" | "LRT" | "p"
'''
def createCSV(fileName):
    file = open(fileName, 'r')
    fileName = fileName.strip(".results")
    newFile = open(f"{fileName}.csv", 'w')
    newFile.write("population,individual,LRT,p\n") # creates header

    for line in file:

        if (not line.startswith("1k")): # keep cycling through lines
            continue
        
        line = line.split() # "1k", "genomes", "group:" "population"
        row  = line[3]

        line = file.readline()
        line = line.split() # "Ancient", "Individual:" "individual"
        row  = f"{row},{line[2]}"

        line = file.readline()
        line = file.readline()
        line = line.strip("\n[]") # "LRT"
        row  = f"{row},{line}"

        line = file.readline()
        line = file.readline()
        line = line.strip("\n[]") # "p value"
        row  = f"{row},{line}\n"

        newFile.write(row)

    file.close()
    newFile.close()

'''
Cleans .results file from running continuity.py.
Removes "Reading line:" parts due to not
being able to write output to a file within
the script. Use this instead of formatTable
if you want to inspect the raw output of
Schraiber's software.
'''
def cleanResults(fileName):
    file = open(fileName, 'r')

    newFile = open(f"cleaned_{fileName}", 'w')

    for line in file:

        if(not line.startswith("Reading line")):
            newFile.write(line)

#############################
#
# Main
#
#############################

fileName = "CHB_CHS_CDX_JPT.results"
createCSV(fileName)