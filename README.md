# DelMissingSite

## 1. Introduction


**DelMissingSite** is designed to automatically remove alignment sites with high proportions of missing data. It can:
- **Process large datasets in parallel**, making it suitable for high-throughput sequencing.
- **Efficiently handle long sequence alignments**, which traditional tools like Gblocks struggle with.

---

## 2. Usage

### 2.1 System Requirements and Dependencies

This software is written in Python 3 and should run on most operating systems, including Linux, Windows, and macOS. It has been successfully tested on Linux and Windows systems. The following dependencies are required:

```bash
Python3
Numpy
Pandas
Biopython
```

### 2.2 Parameters

DelMissingSite is designed to be easy to use, with minimal input required. It automatically reads all `.fasta` files in the current directory as input and only requires two parameters:

```bash
-p  # A floating-point value between 0 and 1 representing the maximum allowed proportion of missing data at a site. Default: 0.2. Sites with missing data exceeding this proportion will be removed.
-n  # An integer (≥1) representing the maximum number of threads to use for parallel processing. Default: 12. The software will process up to 12 alignment files simultaneously.
```

### 2.3 Examples

```bash
python delmissingsite.py -h
# Function: Display the help message.

python delmissingsite.py
# Function: Run the software with default parameters. All `.fasta` files in the current directory will be processed, and any site with >20% missing data will be removed. The software will process up to 12 files in parallel.

python delmissingsite.py -p 0.12 -n 40
# Function: Run the software with a custom threshold for missing data and parallel processing. Any site with >12% missing data will be removed. The software will process up to 40 files in parallel.
```

---

## Citation

If you use this software in your research, please cite the following:

He, J. et al. (2022). A phylotranscriptome study using silica gel-dried leaf tissues produces an updated robust phylogeny of Ranunculaceae. Molecular Phylogenetics and Evolution, 174, 107545.

