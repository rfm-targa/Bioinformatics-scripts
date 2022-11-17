# chewBBACA


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
