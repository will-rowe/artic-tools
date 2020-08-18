# Commands

***

## align_trim

The `align_trim` command is used to softmask primers from an alignment. It mirrors the latest version of the artic pipeline `align_trim` command [found here](https://github.com/artic-network/fieldbioinformatics/blob/master/artic/align_trim.py).

Example usage:

```
artic-tools align_trim -b in.bam primerscheme.bed > out.bam 2> out.log
```

## validate_scheme

The `validate_scheme` command is used to check your primer scheme confirms to an ARTIC standard and can be used in our pipelines.

Example usage:

```
artic-tools validate_scheme primerscheme.bed
```

## check_vcf

The `check_vcf` command is used to check a VCF file and filter variants into PASS and FAIL VCF files.

Example usage:

```
artic-tools check_vcf primerscheme.bed in.vcf
```