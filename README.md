# Bioinformatics-scripts

Collection of Python scripts to perform simple bioinformatics tasks including post-processing of results data from well-known tools.

## Features

The scripts in this repository are organized into the following categories:

- **FASTA**: FASTA file processing for DNA/protein sequence translation, sequence summary statistics, protein diversity analysis.
- **Download**: Download of sequence data from public repositories.
- **SNP**: SNP extraction and analysis.
- **MLST**: Post-processing of Multilocus Sequence Typing (MLST) results.
- **BAM**: BAM file processing utilities.
- **Figures**: Scripts used to create visualizations.

## Installation

The scritps in this repository use the following Python 3.x packages:

- [NumPy](https://numpy.org/)
- [pandas](https://pandas.pydata.org/)
- [Plotly](https://plotly.com/python/)
- [SciPy](https://scipy.org/)
- [Biopython](https://biopython.org/)
- [pyGenomeViz](https://github.com/moshi4/pyGenomeViz)

You can use conda to create an environment with all these packages using the following commands:

```bash
# Create a new environment named bioscripts with Python 3.x
conda create -n bioscripts python
# Activate the environment
conda activate bioscripts
# Install the packages
conda install -c bioconda -c conda-forge numpy pandas plotly scipy biopython pygenomeviz
```
