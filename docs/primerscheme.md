# Primer Scheme

Supported primer schemes are found in the [ARTIC repos](https://github.com/artic-network).

This doc page is a work in progress...

## Primer processing

The following tags are required to exist at the end of the primer IDs:

| tag    | meaning              |
| ------ | -------------------- |
| _LEFT  | the left primer      |
| _RIGHT | the right primer     |
| _alt   | the primer is an alt |

For example:

```
some_id_string_42_LEFT
some_id_string_42_RIGHT
some_id_string_42_RIGHT_alt
```

In the above example, the 3 primers are grouped as one amplicon in the scheme and the two right hand primers will be merged.

To merge a primer, the primers are combined such that a maximal span is achieved. The merged primer will have `_alt` dropped from the ID.

A canonical primer ID is an ID where all tags have been removed. So in the above example, the canonical ID will be `some_id_string_42`.

TODO: more info on primer processing logic


## Scheme Validation

The primer schemes are read from file and validated.

On reading from file, the following must be true:

* file must exist and be readable
* a recognised primer scheme version must be provided
* must be TSV with correct column count correct for scheme version
* must not contain multiple reference sequence IDs
* each row must encode a primer (problem rows are flagged and validation fails after all rows are tried)

Once the file processed for primers, the following checks are made:

* check there are primers in the scheme
* check the number of forward primers match the number of reverse primers (this is **after** merging alts)
* check forward and reverse primers make proper amplicons (based on shared canonical primer IDs)
* check there are no gaps in the scheme

The following information can be reported from the scheme:

* number of pools, primers, amplicons etc.
* mean amplicon size
* the scheme span, with respect to the reference co-ordinates
* the scheme amplicon overlaps (i.e. the proportion of the scheme span with >1 amplicon coverage)
