# Primer Schemes

## About

Both `artic tools` and the `artic pipeline` are designed to work with tiling amplicon schemes for targeted enrichment of viral genomes. The protocol is described in the [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066) Nature Protocols paper and the schemes are produced by [Primal Scheme](https://primalscheme.com/).

Tiling amplicon schemes contain primers that are designed using a greedy algorithm against multiple reference genomes. Tiling refers to the fact that neighbouring amplicons overlap one another, ensuring that complete coverage of the target genome is acheived. To prevent short overlap products being produced between neighbouring amplicons, 2 primer pools are used to alternate primer pairs. The following figure from the [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066) paper shows how primers from 2 pools are used to tile across a reference genome:

![primal-scheme-fig](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnprot.2017.066/MediaObjects/41596_2017_Article_BFnprot2017066_Fig3_HTML.jpg?as=webp)

> [Quick J et al. 2017](https://www.nature.com/articles/nprot.2017.066)
> (a) Schematic showing the regions amplified in pools 1 (upper track) and 2 (lower track), and the intended overlap between pools (as determined in Step 1)...

## Availability

Frequently used primer schemes for viral genome sequencing are found in the [ARTIC primer scheme repo](https://github.com/artic-network/primer-schemes). These include Nipah, Ebola and SARS-CoV-2 primer schemes. Multiple scheme versions may be available, with a higher version number denoting schemes that have had improvements made to them (e.g. introduction of alternate primers).

ARTIC primer schemes can be downloaded using the [artic-tools get_scheme](./commands.md#get_scheme) command. For example:

```sh
artic-tools get_scheme ebola --schemeVersion 2
```

The `artic pipeline` will also attempt to download a scheme using `artic tools` if a scheme can not be found locally.

## File format

Primer schemes produced by [Primal Scheme](https://primalscheme.com/) come with several files (which you can read more about [here](https://github.com/aresti/primalscheme#output)). The 2 files which are needed for `artic tools` are:

| file extension      | description                                                                    |
| ------------------- | ------------------------------------------------------------------------------ |
| `*.reference.fasta` | the sequence of the reference genome used in the scheme                        |
| `*.primer.bed`      | the coordinates of each primer in the scheme, relative to the reference genome |

**Note**: `*.primer.bed` was formerly `*.scheme.bed`, but has lately changed to reflect a slightly updated file format. The old `*.scheme.bed` should still work fine and is available alongside the new versions at the [ARTIC primer scheme repo](https://github.com/artic-network/primer-schemes).

The `*.primer.bed` file is in 6-column BED format, with the following column descriptions:

| column | name       | type   | description                                               |
| ------ | ---------- | ------ | --------------------------------------------------------- |
| 1      | chrom      | string | primer reference sequence                                 |
| 2      | chromStart | int    | starting position of the primer in the reference sequence |
| 3      | chomEnd    | int    | ending position of the primer in the reference sequence   |
| 4      | name       | string | primer name                                               |
| 5      | primerPool | int    | primer pool<sup>\*</sup>                                  |
| 6      | strand     | string | primer direction (+/-)                                    |

<sup>\*</sup> column 5 in the BED spec is an int for score, whereas here we are using it to denote primerPool.

**Note**: BED format, which is a 0-based, half-open format. This means that reference sequence position counting starts at 0 and the chromEnd is not included in the primer sequence.

## Scheme logic

The following sections will discuss how primer schemes are read from a `*.primer.bed` file and used in the ARTIC software.

### Reading schemes

Each line in a `*.primer.bed` file is a single primer. Lines are processed one at a time and the input file does not to be sorted in any way beforehand (although most schemes are sorted by primer start coordinate). As lines are read, they are converted to primer objects and added to an in-memory scheme.

Most of the logic behind creating a primer object relies on the primer name (column 4). As a primer line is read, we check for the following tags in the name field:

| tag      | meaning              |
| -------- | -------------------- |
| `_LEFT`  | the left primer      |
| `_RIGHT` | the right primer     |
| `_alt`   | the primer is an alt |

A `_LEFT` or `_RIGHT` tag is required and only one is allowed per primer. The `_alt` tag is optional and denotes that the primer is an alternate belonging to a primer pair.

**Note**: the tags are cases sensitive.

Depending on the presence of the `_LEFT` or `_RIGHT` tag, the primer is assigned a direction (+/- to denote forward/reverse). This information is also encoded in column 6 of newer schemes. Once the tags have been checked for they are removed and we are left with a primer base ID which is used to group primers by amplicon. For example:

```
REGION_42_LEFT
REGION_42_RIGHT
REGION_42_RIGHT_alt
```

> These primers will be treated as one amplicon in the scheme, which is given the base ID "REGION_42". The two right primers will be merged to yield a single right primer.

To merge a primer with its alternate, the primers are combined such that a maximal span is achieved. The merged primer will have `_alt` dropped from the ID.

In addition to using the name field for creating the primer object, we record the start and end coordinates (columns 2 and 3) and the primer pool (column 5). The primer pool is encoded as an integer in the later scheme versions, but we treat it as a string in the primer scheme code (for backwards compatibility with the older primer scheme files that stored it as a string).

### Scheme validation

The primer schemes are read from file and validated.

On reading from file, the following must be true:

- file must exist and be readable
- a recognised primer scheme version must be provided
- must be TSV with correct column count correct for scheme version
- must not contain multiple reference sequence IDs
- each row must encode a primer (problem rows are flagged and validation fails after all rows are tried)

Once the file processed for primers, the following checks are made:

- check there are primers in the scheme
- check the number of forward primers match the number of reverse primers (this is **after** merging alts)
- check forward and reverse primers make proper amplicons (based on shared canonical primer IDs)
- check there are no gaps in the scheme

The following information can be reported from the scheme:

- number of pools, primers, amplicons etc.
- mean amplicon size
- the scheme span, with respect to the reference co-ordinates
- the scheme amplicon overlaps (i.e. the proportion of the scheme span with >1 amplicon coverage)
