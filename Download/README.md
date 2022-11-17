# Download

### download_ncbi_assemblies.py

#### Description

Accepts a Genome Assembly and Annotation report table from the NCBI and downloads the genome assemblies of the samples listed
in the table. The default behavior will download files from RefSeq. If the genome assembly or file that has to be downloaded
cannot be found in RefSeq, it will try to download from Genbank. The full list of species can be viewed at
https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/ and the script accepts a table with a list of assemblies for a single
species (e.g.: https://www.ncbi.nlm.nih.gov/genome/browse/#!/prokaryotes/186/). Files are download in GZIP format. It is
possible to specify a file extension/suffix to download genome assembly data in a specific file format. You can download the
following file types from the FTP server:

- "assembly_report.txt"
- "assembly_stats.txt"
- "cds_from_genomic.fna.gz"
- "feature_count.txt.gz"
- "feature_table.txt.gz"
- "genomic.fna.gz"
- "genomic.gbff.gz"
- "genomic.gff.gz"
- "genomic.gtf.gz"
- "protein.faa.gz"
- "protein.gpff.gz"
- "rna_from_genomic.fna.gz"
- "translated_cds.faa.gz"

#### Usage

Download genome assemblies in fasta.gz format:

```
python download_ncbi_assemblies.py -t <input_table> -o <output_directory>
```

Download compressed GFF files:

```
python download_ncbi_assemblies.py -t <input_table> -o <output_directory> --fe "genomic.fna.gz"
```

Download only the genome assemblies that are available in RefSeq:

```
python download_ncbi_assemblies.py -t <input_table> -o <output_directory> --ftp refseq
```
