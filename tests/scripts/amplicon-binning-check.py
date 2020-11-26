#!/usr/bin/env python3
import sys
import os

# open up the truth set
dirname = os.path.dirname(__file__)
lookup = {}
with open(os.path.join(dirname, "../data/CVR1.artic.alignreport.txt")) as fh:
    next(fh)
    for line in fh:
        read = line.split("\t")[0]
        assignment = line.split("\t")[3]
        lookup[read] = assignment
fh.close()
numReads = len(lookup)
print("---")
print("reference:\t{}" .format(sys.argv[1]))
print("read total:\t{}" .format(numReads))
print("---")

# open the amplicon binnings
seenReads = {}
correctBinnedReads = {}
incorrectBinnedReads = {}
with open(sys.argv[1]) as fh:
    for line in fh:
        readID = line.split("\t")[0]
        assignment = line.split("\t")[1].rstrip()
        if readID not in lookup:
            print("oh noes, read not found in truth set: {}" .format(readID))
            sys.exit(1)
        seenReads[readID] = 1

        # work out if correct assignment
        if assignment == lookup[readID]:
            correctBinnedReads[readID] = 1
        else:
            # print("{}\tgot:{}\twanted:{}\tkmerFrac:{}" .format(readID, assignment, lookup[readID], line.split("\t")[2].strip()))
            if readID not in incorrectBinnedReads:
                incorrectBinnedReads[readID] = 1
            else:
                incorrectBinnedReads[readID] += 1
fh.close()

# get number of reads which had correct binning AND incorrect binning
intersection = correctBinnedReads.keys() & incorrectBinnedReads.keys()

# create the counts
numCorrectBin = len(correctBinnedReads)
numIncorrectBin = len(incorrectBinnedReads) - len(intersection)
percCorrectBin = (numCorrectBin / numReads) * 100
percIncorrectBin = (numIncorrectBin / numReads) * 100

print("binning file:\t{}" .format(sys.argv[1]))
print("correctly binned:\t{}\t{:.2f}% of reference reads" .format(
    numCorrectBin, percCorrectBin))
print("incorrectly binned:\t{}\t{:.2f}% of reference reads" .format(
    numIncorrectBin, percIncorrectBin))
