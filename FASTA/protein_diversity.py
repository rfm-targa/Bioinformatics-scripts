#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
Given an input directory with a set of FASTA files, this script
translates DNA sequences of all alleles in each input file and
determines which sequences code for the excat same protein
sequence.

This script writes the translated sequences to FASTA files with
sequences ordered by descending length and with a new dot notation
before the allele identifier. The first number is the protein
identifier for which that allele codes and the second number is
the number of the allele that codes for that protein. The notation
'2.3' indicates that the allele codes for the second protein listed
for the gene and that it is the third allele that codes for that
protein.

>seq_id*2 (protein sequence header, the '*' separates the gene
           identifier and the integer identifier given to the
           protein sequence to highlight that it was the second
           distinct protein listed for that gene.)

>seq_id*2.3_9 (dna/allele sequence header, the '*' separates the
               gene identifier from the dot notation and the header
               ends with the original allele identifier.)
"""


import os
import argparse
import multiprocessing

from Bio import SeqIO, Seq
import plotly.graph_objs as go
from plotly.offline import plot


def group_by_protein(fasta_file):
    """Group DNA sequences based on the protein they code for.

    Parameters
    ----------
    fasta_file : str
        Path to the FASTA file with DNA sequences for a gene.

    Returns
    -------
    protein_diversity : dict
        Dictionary with a gene identifier as key and another
        dictionary as value. The nested dictionary has protein
        sequences as keys and a list as value for each key.
        Each list has the allele identifiers and sequences
        that code for that protein, organized in tuples.
    """
    protein_diversity = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seqid = record.id
        allele_id = seqid.split('_')[-1]
        sequence = str(record.seq)
        try:
            protein = Seq.translate(sequence, table=11, cds=True)
        except Exception:
            continue

        protein_diversity.setdefault(protein, []).append((allele_id, sequence))

    return protein_diversity


def attribute_ids(gene_proteins, output_dir, gene_id):
    """Create records for DNA and distinct protein sequences.

    Parameters
    ----------
    gene_proteins : dict
        Dictionary with a gene identifier as key and
        another dictionary as value. The nested
        dictionary has protein sequences as keys and
        a list as value for each key. Each list has
        the allele identifiers and sequences that code
        for that protein, organized in tuples.
    output_dir : str
        Path to the output directory where new files
        will be stored.

    Returns
    -------
    A list with the following variables:
        - the gene identifier;
        - the number of distinct proteins for that gene;
        - the number of distinct DNA sequences for that gene.
    """
    # get proteins and sort by descending length
    gene_prots = gene_proteins
    sorted_proteins = sorted(list(gene_prots.keys()), key=len,
                             reverse=True)

    protid = 1
    alleleid = 1
    protein_counts = {}
    protein_mapping = {}
    dna_records = []
    protein_records = []
    for protein in sorted_proteins:
        current_protein = gene_prots[protein]

        alleles = [a[0] for a in current_protein]
        total_alleles = len(alleles)
        protein_counts[protid] = total_alleles
        protein_mapping[protid] = alleles
        sequences = [a[1] for a in current_protein]

        protein_record = '>{0}*{1}\n{2}'.format(gene_id, protid, protein)
        protein_records.append(protein_record)

        for i, allele in enumerate(alleles):
            dna_sequence = sequences[i]
            dna_record = '>{0}*{1}.{2}_{3}\n{4}'.format(gene_id,
                                                        protid,
                                                        alleleid,
                                                        allele,
                                                        dna_sequence)

            dna_records.append(dna_record)
            alleleid += 1

        alleleid = 1
        protid += 1

    return [protein_mapping, protein_counts,
            protein_records, dna_records]


def write_lines(lines, output_file):
    """Write a list of lines to a file.

    Parameters
    ----------
    lines : list
        List with strings/lines.

    output_file : str
        Path to the output file.
    """
    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(lines)+'\n')

    return output_file


def protein_diversity(file, output_dir):
    """Evaluate protein diversity based on input FASTA file.

    Parameters
    ----------
    file : str
        Path to FASTA file with DNA sequences.
    output_dir : str
        Output directory that will be used to store the
        new files.

    Returns
    -------
    Writes a FASTA file with the input DNA sequences,
    including dot notation in the header and a FASTA file
    with the distinct protein sequences resulting from the
    translation of all DNA sequences, ordered by descending
    length.
    """
    gene_id = os.path.basename(file).split('.fasta')[0]

    try:
        gene_proteins = group_by_protein(file)
        locus_data = attribute_ids(gene_proteins, output_dir, gene_id)
        # write files with data
        dna_file = '{0}/{1}_dna.fasta'.format(output_dir, gene_id)
        write_lines(locus_data[3], dna_file)

        protein_file = '{0}/{1}_protein.fasta'.format(output_dir, gene_id)
        write_lines(locus_data[2], protein_file)

        # write file with protein counts
        protein_counts_file = '{0}/{1}_protein_counts.tsv'.format(output_dir, gene_id)
        with open(protein_counts_file, 'a') as outfile:
            lines = [(k, v) for k, v in locus_data[1].items()]
            lines = sorted(lines, key= lambda x: x[1], reverse=True)
            lines = ['protein_id\tallele_count'] + ['{0}\t{1}'.format(*l) for l in lines]
            lines = '\n'.join(lines)
            outfile.write(lines+'\n')

        print('Processed file: {0}'.format(file))
        return [gene_id, len(locus_data[2]), len(locus_data[3])]
    except Exception:
        print('Failed for file: {0}'.format(file))


def plot_diversity(data, title, output_dir):
    """Create Bar plot for the number of distinct proteins and DNA sequences.

    Parameters
    ----------
    data : list
        List with sublists. Each sublist has
        3 variables: gene identifier, number of distinct
        proteins and number of distinct alleles.
    title : str
        Title for the Bar plot.
    output_dir : str
        Output diretory where the HTML file will
        be created.

    Returns
    -------
    Creates a HTML file that can be opened with
    a Browser to visualize a plot comparing the
    number of distinct alleles and the number of
    distinct proteins for each gene.
    """
    y_alleles = [d[2] for d in data]
    y_proteins = [d[1] for d in data]
    x_data = [d[0] for d in data]

    alleles_tracer = go.Bar(name='#alleles', x=x_data,
                            y=y_alleles, marker_color='rgb(56,108,176)')
    protein_tracer = go.Bar(name='#proteins', x=x_data,
                            y=y_proteins, marker_color='rgb(127,201,127)')

    fig = go.Figure(data=[alleles_tracer, protein_tracer])

    complete_title = '{0}<br>Distinct alleles vs distinct proteins'.format(title)
    fig.update_layout(barmode='group', title_text=complete_title)

    fig.update_xaxes(title_text='Gene identifier')
    fig.update_yaxes(title_text='Number of distinct alleles<br>& '
                    'Number of distinct proteins')

    output_file = os.path.join(output_dir, '{0}.html'.format(title))

    plot(fig, filename=output_file, auto_open=False)


def main(input_directory, output_directory, threads, plot_creation):

    # get list of schema files
    schema_files = [os.path.join(input_directory, file)
                    for file in os.listdir(input_directory)
                    if '.fasta' in file]
    schema_basename = os.path.basename(input_directory) \
                      if input_directory[-1] != '/' \
                      else os.path.basename(input_directory[:-1])

    # create directory to store new files
    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    # cluster alleles identifiers based on protein sequence
    plot_data = []
    pool = multiprocessing.Pool(threads)
    for file in schema_files:
        pool.apply_async(protein_diversity,
                         (file, output_directory),
                         callback=plot_data.append)
    pool.close()
    pool.join()

    if plot_creation is True:
        plot_diversity(plot_data, schema_basename, output_directory)

    print('Done!')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-directory', type=str,
                        required=True, dest='input_directory',
                        help='The chewBBACA schema directory.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Output directory that will be created '
                             'to store the new files.')

    parser.add_argument('-t', '--threads', type=int,
                        required=False, default=1, dest='threads',
                        help='Number of CPU cores to use.')

    parser.add_argument('-p', '--plot-creation', action='store_true',
                        dest='plot_creation',
                        help='If the process should create a Bar plot '
                             'comparing the number of alleles and the '
                             'number of distinct proteins for each gene.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
