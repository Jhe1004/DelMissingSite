# DelMissingSite

## 1. Introduction

### 1.1 Purpose

The diversity of life forms on Earth is a reflection of their shared evolutionary history. Understanding the relationships between different biological groups is not only essential for evolutionary biology but also forms the foundation for taxonomy, species classification, and various biological disciplines. A key focus in these studies is the construction of reliable phylogenetic trees, which represent the evolutionary relationships between species.

With advancements in molecular biology, researchers now rely on macromolecules like DNA and protein sequences to infer evolutionary histories more accurately (Yang & Olmstead, 1997; Qiu et al., 2000; Chase et al., 2016). The typical workflow of phylogenetic research involves two major steps: **sequence alignment** and **phylogenetic tree construction**. Various tools have been developed to aid in these processes (Thompson et al., 1994; Wilgenbusch & Swofford, 2003; Edgar et al., 2004; Darling et al., 2004; Ronquist et al., 2012; Katoh et al., 2013; Stamatakis, 2014; Minh et al., 2020).

However, some alignment sites may not be strictly orthologous due to significant genetic distance between sequences, which could impact the accuracy of the phylogenetic tree. Such sites are often characterized by excessive gaps (Löytynoja & Goldman, 2008). Therefore, a crucial yet underexplored step between sequence alignment and tree construction is the **removal of sites with excessive missing data (gaps)** based on a user-defined threshold.

Historically, due to shorter sequence lengths, researchers manually removed these gap-rich sites. Automation tools like Gblocks (Castresana, 2000) were later developed to handle this process. However, with the advent of high-throughput sequencing technologies, researchers now deal with thousands of sequences or even whole genomes, making manual deletion impossible. Gblocks also has limitations: it is single-threaded and struggles with long sequences, leading to significant time delays. Thus, there is a pressing need for a tool that can efficiently handle large-scale and long-sequence datasets in parallel, similar to Gblocks, but with higher scalability and efficiency.

### 1.2 Features

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

He J. (2022). Jhe1004/DelMissingSite: (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.6415293

---

### Summary of Enhancements in This Version:
- **Parallel Processing**: DelMissingSite can handle multiple alignment files simultaneously, making it suitable for large-scale datasets.
- **Efficiency**: The software can process long sequences much faster than traditional tools.
