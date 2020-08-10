# Softmasker

The artic pipeline softmasker is used to trim alignments within amplicon boundaries, as well as perform some basic filtering and normalisation. The original softmasker code can be found in the field bioinformatics repository, within the `align_trim` module [here](https://github.com/artic-network/fieldbioinformatics/blob/master/artic/align_trim.py).

## Original softmasker

The basic workflow of the original softmasker is as follows:

### input

* an alignment of amplicon sequencing data to a reference genome
* a primer scheme (in ARTIC BED format)

### output

* a filtered alignment of trimmed alignment segments
* a logfile/report

### process

* filter each alignment segment from an input alignment to remove unmapped or supplementary segments
* locate the nearest forward and reverse primers for each segment
* add a read group to the segment to denote a primer pool or unmatched
* softmask the segment to be within an amplicon boundary (determined by the primer pair)
* filter the output alignment to remove segments where softmasking has made the alignment ambiguous

### options

* remove alignment segments with incorrectly paired primers
* don't assign read groups with primer pool information
* instead of softmasking segments to include primer sequences, softmask to exclude them
* normalise amplicon coverage by dropping alignment segments after a threshold is reached

## Improved (?) softmasker

The same inputs and workflow are used, but we have achieved some performance improvements (from shuffling logic and porting to CPP) and have made a few minor changes:

* optionally filter segments by mapping quality
* perform both softmasks in single BAM iteration (results in multiple BAM outputs)
* during primer search, in the case that a given position is equidistant between two primer sites the upper bound is now used instead of the lower when locating the primer start sites (this typically ends up with more correctly paired primers)
* output more stats at the end of the report file

The pseudocode for the artic-tools softmasker is:

```
for segment in alignment:

    # filter segment
    if (segment == unmapped) || (segment == supplementary):
        continue
    if (segment.mappingQuality < threshold):
        continue

    # get predicted amplicon for segment 
    amplicon = getAmplicon(segment)
    if (amplicon != properlyPaired):
        continue

    # check abundance
    amplicon.abundance++
    if (amplicon.abundance > threshold):
        continue
    
    # assign primer pool to segment
    assign_readgroup(segment, amplicon)

    # softmask segment within amplicon boundary
    softmask(segment, amplicon, inclPrimers)
    write(outfile1, segment)

    # softmask segment further to exclude primers
    softmask(segment, amplicon, exclPrimers)
    write(outfile2, segment)

``` 

## Not yet implemented in the release candidate

* option to output short fragments in a separate file or in the report
* additional BAM workers
* amplicon coverage and plotting (implemented in primer scheme logic)
* diversity profile for later read group var matching
* overlap region identification (implemented in primer scheme logic)

## Notes

* if normalisation is selected, the threshold is applied in each alignment direction (behaviour the same in the original artic codebase)
* non-ARTIC primer schemes aren't validated against this program
