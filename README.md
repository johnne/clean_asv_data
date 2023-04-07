# Clean ASV data
Code to clean up ASV clustering results and generate stats

## Installation

- Clone the repository and change directory into `clean_asv_data`
```bash
git clone git@github.com:johnne/clean_asv_data.git 
cd clean_asv_data
```

- Create the environment with [conda](https://docs.conda.io/en/latest/miniconda.html)
or [mamba](https://github.com/conda-forge/miniforge#mambaforge):

```bash
mamba env create -f environment.yml
```

- Activate the environment

```bash
conda activate clean_asv_data
```

- Install the scripts with `pip`

```
python -m pip install .
```

You may have to deactivate and re-activate the `clean_asv_data` environment 
after this command.

## Usage

### rename-samples

The `rename-samples` script allows you to rename sample names in a 
tab-separated file where each sample is a column. It uses regular 
expressions defined with the `--regex` argument. Any number of regular 
expressions can be supplied, in the format `<pattern><split><repl>` where 
`<pattern>` is the string pattern to search for, `<repl>` is the string to 
replace the pattern with and `<split>` is a string used to separate the two 
former.

Example:

With a file `counts.tsv´ with sample names in columns like so:

```
        FL001_L1001     FL001_neg       SOIL_1001
ASV1    ...
```

To rename all samples starting with `FL\d+_L` ('FL' followed by 1 or more 
digits (0-9), and the '_L' string) and keep only the `L`:

```bash
rename-samples --regex 'FL\d+_L,L' --regex-split ',' counts.tsv > counts.renamed.tsv
```

Result:
```
        L1001   FL001_neg       SOIL_1001
ASV 
```

All arguments:
```bash
usage: rename-samples [-h] [--regex REGEX [REGEX ...]] [--regex-split REGEX_SPLIT] [--chunksize CHUNKSIZE] input

positional arguments:
  input                 Tab-separated input file with sample names in columns

options:
  -h, --help            show this help message and exit
  --regex REGEX [REGEX ...]
                        One or more regular expressions to use to rename column names of the input file.
  --regex-split REGEX_SPLIT
                        Character used to split the regular expressions into<pattern> and <repl>. For example with --regex 'FL\d+_L,L the --regex-split ',' will replace 'FL\d+_L' with 'L'. Default ','
  --chunksize CHUNKSIZE
                        If input file is very large, specify chunksize to read it in a number of lines at a time

```

### clean-asv-data

The `clean-asv-data` script can be used to clean up an ASV cluster file 
based on several parameters.

Example:

```bash
clean-asv-data --counts asv_counts.tsv --taxonomy asv_taxonomy.tsv 
--clean_rank Family > asv_taxonomy.cleaned.tsv
```

All arguments:
```bash
usage: 
        This script cleans clustering results by removing ASVs if:
        - unassigned or ambiguous taxonomic assignments (e.g.  
        'unclassified' or '_X' in rank labels) 
        - if belonging to clusters present in > max_blank_occurrence% of blanks
        - if belonging to clusters with < min_clust_count total reads
        
       [-h] [--counts COUNTS] [--taxonomy TAXONOMY] [--blanks BLANKS] [--output OUTPUT] [--clean_rank CLEAN_RANK] [--max_blank_occurrence MAX_BLANK_OCCURRENCE] [--blank_removal_mode {cluster,asv}]
       [--min_clust_count MIN_CLUST_COUNT] [--chunksize CHUNKSIZE] [--nrows NROWS]

options:
  -h, --help            show this help message and exit

input/output:
  --counts COUNTS       Counts file of ASVs
  --taxonomy TAXONOMY   Taxonomy file for ASVs. Should also include acolumn with cluster designation.
  --blanks BLANKS       File with samples that are 'blanks'
  --output OUTPUT       Output file with cleaned results

params:
  --clean_rank CLEAN_RANK
                        Remove ASVs unassigned at this taxonomic rank (default Family)
  --max_blank_occurrence MAX_BLANK_OCCURRENCE
                        Remove ASVs occurring in clusters where at least one member is present in <max_blank_occurrence>% of blank samples. (default 5)
  --blank_removal_mode {cluster,asv}
                        How to remove sequences based on occurrence in blanks. If 'asv' (default) remove only ASVs that occur in more than <max_blank_occurrence>% of blanks. If 'cluster', remove ASVs
                        in clusters where one or more ASVs is above the <max_blank_occurrence> threshold
  --min_clust_count MIN_CLUST_COUNT
                        Remove clusters with < <min_clust_count> summed across samples (default 3)

debug:
  --chunksize CHUNKSIZE
                        Size of chunks (in lines) to read from countsfile
  --nrows NROWS         Rows to read from countsfile (for testing purposes only)
```

### generate-statsfile

The `generate-statsfile` script calculates total sum and occurrence of ASVs 
using a countsfile as input.

All arguments:
```bash
usage: generate-statsfile [-h] countsfile

positional arguments:
  countsfile  ASV counts file. Tab-separated, samples in columns, ASVs in rows

options:
  -h, --help  show this help message and exit
```

### count-clusters

The `count-clusters` script sums read counts for each ASV cluster and writes 
to standard out.

All arguments:
```bash
usage: count-clusters [-h] [--clust_column CLUST_COLUMN] [--chunksize CHUNKSIZE] countsfile clustfile

positional arguments:
  countsfile            Tab-separated file with counts of ASVs (rows) in samples (columns)
  clustfile             Tab-separated file with ASV ids in first column and a column specifying the cluster it belongs to

options:
  -h, --help            show this help message and exit
  --clust_column CLUST_COLUMN
                        Name of cluster column. Defaults to 'cluster'
  --chunksize CHUNKSIZE
                        If countsfile is very large, specify chunksize to read it in a number of lines at a time
```

