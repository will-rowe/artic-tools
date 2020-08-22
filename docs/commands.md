# Commands

---

## align_trim

The `align_trim` command is used to softmask primers from an alignment. It mirrors the latest version of the artic pipeline `align_trim` command [found here](https://github.com/artic-network/fieldbioinformatics/blob/master/artic/align_trim.py).

Example usage:

```
artic-tools align_trim -b in.bam primerscheme.bed > out.bam 2> out.log
```

## get_scheme

The `get_scheme` command can download primer schemes and sequences for several ARTIC references.

Example usage:

```
artic-tools get_scheme nipah
```

## validate_scheme

The `validate_scheme` command is used to check your primer scheme confirms to an ARTIC standard and can be used in our pipelines.

Example usage:

```
artic-tools validate_scheme primerscheme.bed
```

It reports some basic stats and can also be used to produce a multifasta of all your primer sequences. Example output looks like this:

```
primer scheme file:     myscheme.bed
primer scheme version:  3
reference sequence ID:  MN908947.3
number of pools:        2
number of primers:      218 (includes 22 alts)
number of amplicons:    98
mean amplicon size:     393
scheme ref. span:       30-29866
scheme overlaps:        29.1326%
primer sequences:       primers.fasta
```

## check_vcf

The `check_vcf` command is used to check a VCF file and to (optionally) filter variants into a PASS VCF file.

Example usage:

```
artic-tools check_vcf --dropPrimerVars --dropOverlapFails -o pass.vcf primerscheme.bed in.vcf
```

Example ouput:

```
[14:13:50] artic-tools::vcfchecker: starting VCF checker
[14:13:50] artic-tools::vcfchecker:     filtering variants: true
[14:13:50] artic-tools::vcfchecker:     output file: pass.vcf
[14:13:50] artic-tools::vcfchecker:     discard primer site vars: true
[14:13:50] artic-tools::vcfchecker:     discard overlap fail vars: true
[14:13:50] artic-tools::vcfchecker: variant at pos 241: C->T
[14:13:50] artic-tools::vcfchecker: variant at pos 3037: C->T
[14:13:50] artic-tools::vcfchecker: variant at pos 12733: C->T
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     nothing seen at position yet, holding var
[14:13:50] artic-tools::vcfchecker: variant at pos 12733: C->T
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     multiple copies of var found at pos 12733 in overlap region
[14:13:50] artic-tools::vcfchecker: variant at pos 14408: C->T
[14:13:50] artic-tools::vcfchecker: variant at pos 22863: TA->T
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     nothing seen at position yet, holding var
[14:13:50] artic-tools::vcfchecker: variant at pos 22868: TG->T
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     var pos does not match with that of previously identified overlap var, holding new var (and dropping held var at 22862)
[14:13:50] artic-tools::vcfchecker: variant at pos 22896: T->TTGG
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     var pos does not match with that of previously identified overlap var, holding new var (and dropping held var at 22867)
[14:13:50] artic-tools::vcfchecker: variant at pos 22909: TA->T
[14:13:50] artic-tools::vcfchecker: variant at pos 22913: T->C
[14:13:50] artic-tools::vcfchecker: variant at pos 22916: CT->TC
[14:13:50] artic-tools::vcfchecker: variant at pos 22926: T->TAA
[14:13:50] artic-tools::vcfchecker: variant at pos 22948: ACC->A
[14:13:50] artic-tools::vcfchecker: variant at pos 22995: C->CA
[14:13:50] artic-tools::vcfchecker: variant at pos 22997: C->T
[14:13:50] artic-tools::vcfchecker: variant at pos 23009: G->GT
[14:13:50] artic-tools::vcfchecker: variant at pos 23057: C->A
[14:13:50] artic-tools::vcfchecker: variant at pos 23098: AC->GT
[14:13:50] artic-tools::vcfchecker: variant at pos 23183: T->TC
[14:13:50] artic-tools::vcfchecker:     located within an amplicon overlap region
[14:13:50] artic-tools::vcfchecker:     var pos does not match with that of previously identified overlap var, holding new var (and dropping held var at 22895)
[14:13:50] artic-tools::vcfchecker: variant at pos 23403: A->G
[14:13:50] artic-tools::vcfchecker: variant at pos 27752: C->T
[14:13:50] artic-tools::vcfchecker: variant at pos 28881: GGG->AAC
[14:13:50] artic-tools::vcfchecker: finished checking
[14:13:50] artic-tools::vcfchecker:     dropped var at 23182 which is in an amplicon overlap region but was only found once
[14:13:50] artic-tools::vcfchecker:     22 variant records processed
[14:13:50] artic-tools::vcfchecker:     18 variant records passed checks
```
