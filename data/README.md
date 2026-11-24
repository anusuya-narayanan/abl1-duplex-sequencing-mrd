# Input data

The analysis script needs the two **Variant Annotator** tables from the Galaxy history.

| File name | Galaxy dataset | What it is |
|---|---|---|
| `dcs_variant_annotator.tsv` | 20, "Variant Annotator on dataset 19" | Per-position counts from duplex consensus reads |
| `sscs_variant_annotator.tsv` | 27, "Variant Annotator on dataset 26" | Per-position counts from single-strand consensus reads |

Each table has one row per chr9 position with columns
`SAMPLE, CHR, POS, +A, +C, +G, +T, -A, -C, -G, -T, CVRG, ALLELES, MAJOR, MINOR, MAF, BIAS`.
`+` and `-` are reads mapped in the forward and reverse direction, `CVRG` is depth and `MAF` is the minor allele frequency.

Raw reads (SRA SRR1799908) and BAM files are not stored here; they are public and large.
