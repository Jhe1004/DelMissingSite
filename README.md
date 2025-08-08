# DelMissingSite

`DelMissingSite` is a high-performance Python script designed to efficiently filter sites with a high proportion of missing data from multiple sequence alignments. This process, often called "alignment trimming" or "filtering," is a crucial step in preparing data for robust phylogenetic analysis.

The script leverages Python's `multiprocessing` module to execute tasks in parallel, significantly speeding up the processing of large numbers of alignment files. It features an automated, multi-stage workflow that can handle extremely long alignments by splitting them into manageable chunks, processing them, and re-concatenating them seamlessly.

---

## Table of Contents
- [Key Features](#key-features)
- [Workflow Overview](#workflow-overview)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [How to Run](#how-to-run)
- [Command-Line Options](#command-line-options)
- [Input and Output](#input-and-output)

---

## Key Features

-   **Site-based Filtering**: Removes alignment columns (sites) where the percentage of missing data (e.g., `-`, `?`, `N`) exceeds a user-defined threshold.
-   **Parallel Processing**: Utilizes multiple CPU cores (`-n` option) to process many files or file chunks simultaneously, drastically reducing computation time.
-   **Automated Large Alignment Handling**: If an input alignment is longer than a built-in threshold (50,000 bp), the script will automatically split it, process the chunks in parallel, and correctly reassemble the filtered result.
-   **Comprehensive Data Cleaning**:
    -   Standardizes multi-line FASTA files into a single-line format.
    -   Removes uninformative sequences (rows) that consist entirely of missing data.
    -   Pads shorter sequences within an alignment to ensure all sequences have uniform length.
-   **Batch Operation**: Automatically detects and processes all `.fasta` files located in the script's current working directory.

## Workflow Overview

The script operates in a three-stage, fully automated pipeline. Understanding this helps interpret the intermediate files you might see during its execution.

1.  **Stage 1: Preprocessing (`.fasta` → `.fa`)**
    -   The script first finds all input files ending in `.fasta`.
    -   It cleans and standardizes each file (format conversion, removal of empty sequences).
    -   If a file is very long, it's split into smaller chunks.
    -   The output of this stage are intermediate files ending in `.fa`.

2.  **Stage 2: Filtering (`.fa` → `.fas`)**
    -   The script takes the `.fa` files and applies the core filtering logic.
    -   It calculates the proportion of missing data for each column and removes columns that exceed the threshold set by the `-p` option.
    -   The filtered results are written to files ending in `.fas`. If the original file was split, there will be multiple `.fas` files for it.

3.  **Stage 3: Concatenation (re-assembly)**
    -   This final stage activates only if files were split in Stage 1.
    -   It identifies all filtered chunks belonging to the same original alignment.
    -   It concatenates these chunks back into a single, complete `.fas` file, ensuring all sequences are correctly ordered and padded.
    -   All temporary intermediate files (`.fa`, `.split.*`, `.temp`) are deleted upon completion.

## Dependencies

The script requires the following Python libraries:
-   `NumPy`
-   `Pandas`
-   `Biopython`

## Installation

Install the required libraries using pip:

    pip install numpy pandas biopython

## How to Run

1.  **Prepare your files**: Place the `DelMissingSite.py` script and all your input `.fasta` alignment files into the same directory.
2.  **Open your terminal**: Navigate to this directory.
3.  **Execute the script**: Run the script from the command line, adjusting the options as needed.

**Example 1: Using default settings**
This will filter all `.fasta` files, removing sites with 20% or more missing data, using 40 CPU cores.

    python DelMissingSite.py

**Example 2: Using custom settings**
This will remove sites with 50% or more missing data, using 16 CPU cores.

    python DelMissingSite.py -p 0.5 -n 16

## Command-Line Options

| Argument | Description |
|---|---|
| `-p`, `--proportion` | The maximum allowed proportion of missing data at any given site. Sites exceeding this proportion will be deleted. **Default: `0.2`**. |
| `-n`, `--num_cpu` | The number of CPU cores to use for parallel processing. Higher values speed up the script when you have many files. **Default: `40`**. |

## Input and Output

-   **Input**: The script automatically finds and processes all files in the current directory that end with the `.fasta` extension.

-   **Output**: For each input file (e.g., `gene1.fasta`), the script will produce a single, filtered output file with the `.fas` extension (e.g., `gene1.fas`). All intermediate files are automatically cleaned up.
