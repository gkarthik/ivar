Manual {#manualpage}
============

## Available Commands

| Command | Description |
|:--------|:------------|
| trim | Trim reads in aligned BAMs |
| variants | Call variants from aligned BAM file |
| filtervariants | Filter variants across replicates
| consensus | Call consensus from aligned BAM file |
| getmasked | Get amplicons with primer mismatches |
| removereads | Remove reads from trimmed BAM file |
| version | Show version information |
| trimadapter | (EXPERIMENTAL) Trim adapter sequences from reads |

To view detailed usage for each command type `ivar <command>`
Note : Commands maked (EXPERIMENTAL) are still under active development.

## Description of commands

#### Trim primer sequences with iVar

Command:
```
ivar trim

Usage: ivar trim -i <input.bam> -b <primers.bed> -p prefix [-m <min-length>] [-q <min-quality>] [-s <sliding-window-width>]

Input Options    Description
           -i    (Required) <input.bam> Indexed aligned bam file to trim primers and quality
           -b    (Required) <primers.bed> BED file with primer sequences and positions
           -m    Minimum length of read to retain after trimming (Default: 30)
           -q    Minimum quality threshold for sliding window to pass (Default: 20)
           -s    Width of sliding window (Default: 4)

Output Options   Description
           -p    (Required) Prefix for the output BAM file
```

Example Usage:
```

```

#### Call variants with iVar

Command:
```
ivar variants

Usage: samtools mpileup -A -d 300000 --reference <reference-fasta> -Q 0 -F 0 <input-bam> | ivar variants -p prefix [-q <min-quality>] [-t <min-frequency-threshold>]

Note : samtools mpileup output must be piped into ivar variants

Input Options    Description
           -q    <min-quality> Minimum quality threshold to count base (Default: 20)
           -t    <min-frequency-threshold> Minimum frequency threshold(0 - 1) to call variants (Default: 0.03)

Output Options   Description
           -p    (Required) Prefix, prefix for the output tsv variant file
```

Example Usage:
```

```

#### Filter variants across replicates with iVar

Command:
```
ivar filtervariants
Usage: ivar filtervariants [-p prefix] replicate-one.tsv replicate-two.csv ...

Output Options   Description
           -p    (Required) Prefix, prefix for the filtered tsv file

```

Example Usage:
```

```

#### Call a consensus sequences from an aligned BAM file.

Command:
```
ivar consensus
Usage: samtools mpileup -A -d 300000 -Q 0 -F 0 [<input-bam>] | ivar consensus [-p prefix]

Note : samtools mpileup output must be piped into "ivar consensus"

Input Options    Description
           -q    Minimum quality threshold to count base (Default: 20)
           -t    Threshold(0 - 100) to call consensus (Default: 0)
Output Options   Description
           -p    (Required) Prefix, prefix for the output tsv file
```

Example Usage:
```

```

#### Get primers with mismatches to the reference sequence.

Command:
```
ivar getmasked
Usage: ivar getmasked [-i <input-filtered-tsv>] [-b <bed-file>]

Input Options    Description
           -i    (Required) Input filtered variants tsv generated from "ivar filtervariants"
           -b    (Required) Bed file with primer indices
```

Example Usage:
```

```

#### Remove reads associated with mismatched primer indices

Command:
```
ivar removereads
Usage: ivar removereads [-i <input-bam>] [-p prefix] primer-index-1 primer-index-2 primer-index-3 primer-index-4 ...

Input Options    Description
           -i    (Required) Input BAM file run through `ivar trim` command whcih adds the primer number to BAM auxillary data
Output Options   Description
           -p    (Required) Prefix, prefix for the filtered BAM file
```

Example Usage:
```

```
