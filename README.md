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

## Typical workflow

For example, with input files under `data/`:

```bash
data
├── asv_counts.tsv
├── blanks.txt
└── clustfile.tsv
```

Create directory to hold processed output:

```bash
mkdir results
```

### Step 1. Clean ASV data

```bash
clean-asv-data --countsfile data/asv_counts.tsv \
  --blanksfile data/blanks.txt \
  --clustfile data/clustfile.tsv > results/cleaned_clustfile.tsv
```

### Step 2. Generate stats

```bash
generate-statsfile --countsfile data/asv_counts.tsv \
  --blanksfile data/banks.txt > resultsd/asv_stats.tsv
```

### Step 3. Generate consensus taxonomy

```bash
consensus-taxonomy --countsfile data/asv_counts.tsv \
  --clustfile results/cleaned_clustfile.tsv > results/cleaned_cluster_taxonomy.tsv
```

### Step 4. Count clusters

```bash
count-clusters --countsfile data/asv_counts.tsv \
  --clustfile results/cleaned_clustfile.tsv > results/cleaned_cluster_count.tsv
```

## Scripts

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

## consensus-taxonomy

All arguments:
```
usage: consensus-taxonomy [-h] [--countsfile COUNTSFILE] [--clustfile CLUSTFILE] [--configfile CONFIGFILE] [--ranks RANKS [RANKS ...]] [--clust_column CLUST_COLUMN]
                          [--consensus_threshold CONSENSUS_THRESHOLD] [--consensus_ranks CONSENSUS_RANKS [CONSENSUS_RANKS ...]] [--chunksize CHUNKSIZE]

options:
  -h, --help            show this help message and exit
  --countsfile COUNTSFILE
                        Counts file of ASVs
  --clustfile CLUSTFILE
                        Taxonomy file for ASVs. Should also include a column with cluster designation.
  --configfile CONFIGFILE
                        Path to a yaml-format configuration file. Can be used to set arguments.
  --ranks RANKS [RANKS ...]
                        Ranks to include in the output.
  --clust_column CLUST_COLUMN
                        Name of cluster column, e.g. 'cluster'
  --consensus_threshold CONSENSUS_THRESHOLD
                        Threshold (in %) at which to assign taxonomy to a cluster
  --consensus_ranks CONSENSUS_RANKS [CONSENSUS_RANKS ...]
                        Ranks to use for calculating consensus. Must be present in the clustfile.
  --chunksize CHUNKSIZE
                        If countsfile is very large, specify chunksize to read it in a number of lines at a time
```