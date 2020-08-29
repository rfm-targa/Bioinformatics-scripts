# BioinfUtils

Collection of Python scripts to perform various Bioinformatics tasks.

### chewie_allelecall_simdiff.py

#### Description

Accepts a matrix with results from the AlleleCall process of [chewBBACA](https://github.com/B-UMMI/chewBBACA) and 
can apply two processes to that matrix: 

  - mask certain matrix elements by substituting those elements with other characters; 
  - determine the number of shared loci and the number of allele differences between each pair of samples represented
    in the AlleleCall matrix (with the option to mask matrix elements before determining the number of shared alleles 
    and differences).

The default masking option will substitute all ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5, LOTSC and LNF cases with '0' and
the 'INF-' prefix of inferred alleles will always be removed to homogenize valid allele identifiers. Passing a single
word will change the default substitution value for that word.To change specific matrix elements, the string should
be formatted as:

  'ASM=short,ALM=large,LNF=pop'

which will change the default substitution value for ASM and ALM cases, maintaining the default value of '0' to substitute
remaining cases. The 'pop' expression serves to signal that a type of case should not be substituted. 

The ',' character is used to separate several substitution assignments and should only be used that way. The '=' character
is used to change the default substitution value for the right side matrix element (ASM=short --> ASM cases will be substituted
with the 'short' word) and should only be used for that purpose.

#### Usage

- 'mask' process:

```
python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> -p mask -o <masked_output_matrix> -mc ASM=short,ALM=large,LNF=pop
```

Will substitute ASM elements with 'short' and ALM elements with 'large'. LNF cases will not be substituted due to the use
of the 'pop' expression. NIPH, NIPHEM, PLOT3, PLOT5 and LOTSC cases will be substituted by the default '0' character. All
'INF-' prefixes will be removed from inferred alleles.

```
python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> -p mask -o <masked_output_matrix> -mc sub
```

Will change the default subtitution value from '0' to 'sub' and all ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5, LOTSC and LNF
cases will be substituted by 'sub'.

- 'sims' process:

```
python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> -p sims -o <masked_output_matrix> -m True -mc ASM=short,LNF=pop
```

Will substitute ASM elements with 'short'. LNF cases will not be substituted due to the use of the 'pop' expression. ALM,
NIPH, NIPHEM, PLOT3, PLOT5 and LOTSC cases will be substituted by the default '0' character. All 'INF-' prefixes will be
removed from inferred alleles. After masking matrix values, the script will determine a new matrix that has the number of
loci shared between each pair of samples above the diagonal and the number of allele differences below the diagonal. Masked
values will not be considered when determining similarity, only cases where each sample has a valid allele identifier will
be considered.

```
python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> -p sims -o <masked_output_matrix> -m False
```

Will not mask any matrix values due to the 'False' value passed to the 'm' argument. The number of sample similarities and differences will 
be determined based on the full set of matrix elements because no matrix value will be considered as masked.

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

### fasta_integer_ids.py

#### Description

Reads lines from a FASTA file given as input and substitutes the identifiers of the headers by a single integer, starting
at the user defined value (default=1) and incrementing that value by a gap value (default=1) to assign a different integer
to the next sequence.

#### Usage

Default parameters, starting at 1 and incrementing by 1 (1, 2, 3, 4, ...):

```
python fasta_integer_ids.py -i <input_fasta> -o <output_fasta>
```

Start at 3 and incrementing by 2 (3, 5, 7, 9, ...):

```
python fasta_integer_ids.py -i <input_fasta> -o <output_fasta> -si 3 -ig 2
```

### generate_snippy_aln_heatmap.py

#### Description

Creates a basic Heatmap from a [Snippy](https://github.com/tseemann/snippy) core.aln or core.full.aln file.
    
Each sample nucleotide position is colored according to the same position in the reference sequence.

Colorscale hex codes:  
Equal to reference: `#ece7f2` ![#ece7f2](https://via.placeholder.com/15/ece7f2/000000?text=+)  
SNP: `#41ab5d` ![#41ab5d](https://via.placeholder.com/15/41ab5d/000000?text=+)  
Zero coverage or deletions: `#0868ac` ![#0868ac](https://via.placeholder.com/15/0868ac/000000?text=+)  
Low coverage, masked region on reference and heterozygous or poor quality genotype: `#252525` ![#252525](https://via.placeholder.com/15/252525/000000?text=+) 

This script requires that 'Biopython' and 'Plotly' be installed.

Note: the output file size will vary according to the number of samples. When using a core.full.aln file, the output
HTML file might be large depending on the number of samples (e.g: a core.full.aln file with the alignment between the
assembly of 1 reference and 7 samples of Streptococcus agalactiae originated a HTML file with 124MB).

#### Usage

Create basic heatmap based on a core.aln or core.full.aln file generated by Snippy:

```
python generate_snippy_aln_heatmap.py -i <snippy_aln_file> -rid <reference_id> -pf <output_html>
```

