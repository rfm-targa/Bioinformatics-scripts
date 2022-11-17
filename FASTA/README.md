# FASTA


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
