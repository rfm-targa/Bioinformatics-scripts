#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

    Given an input schema, and for each gene file, translates DNA
    sequences of all alleles in each locus and determines 
    which alleles code for the excat same protein sequence.
    
    This script outputs new schema files with sequences ordered by 
    descending length and with a new dot notation before the allele
    identifier. The first algarism is the protein identifier for which
    that allele codes and the second algarism is the number of the 
    allele that codes for that protein. The notation '2.3' indicates
    that the allele codes for the second protein listed for the gene
    and that it is the third allele that codes for that protein.
    
    >gene_id*2 (protein sequence header, the '*' separates the gene
                identifier and the integer identifier given to the
                protein sequence to highlight that it was the second
                distinct protein listed for that gene.)
    
    >gene_id*2.3_9 (dna/allele sequence header, the '*' separates the
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
    """ Groups DNA sequences based on the protein they code for.

        Args:
            fasta_file (str): path to the FASTA file with DNA
            sequences for a gene.

        Returns:
            protein_diversity (dict): dictionary with a gene
            identifier as key and another dictionary as value.
            The nested dictionary has protein sequences as keys
            and a list as value for each key. Each list has
            the allele identifiers and sequences that code for
            that protein, organized in tuples.
    """

    protein_diversity = {}
    basename = os.path.basename(fasta_file)
    protein_diversity[basename] = {}
    for record in SeqIO.parse(fasta_file, 'fasta'):
        seqid = record.id
        allele_id = seqid.split('_')[-1]
        sequence = str(record.seq)
        protein = Seq.translate(sequence, table=11, cds=True)

        if protein in protein_diversity[basename]:
            protein_diversity[basename][protein][0].append((allele_id, sequence))
        else:
            protein_diversity[basename][protein] = [[(allele_id, sequence)]]

    return protein_diversity


def write_files(gene_proteins, output_dir):
    """ Writes a FASTA file with the allele sequences for a gene,
        including dot notation in the sequence headers. Also writes
        a FASTA file with all distinct proteins for that gene.

        Args:
            gene_proteins (dict): dictionary with a gene
            identifier as key and another dictionary as value.
            The nested dictionary has protein sequences as keys
            and a list as value for each key. Each list has
            the allele identifiers and sequences that code for
            that protein, organized in tuples.
            output_dir (str): path to the output directory where
            new files will be stored.

        Returns:
            A list with the following variables:
                - the gene identifier;
                - the number of distinct proteins for that gene;
                - the number of distinct DNA sequences for that gene.
    """

    for g, v in gene_proteins.items():
        gene_id = g.split('_')[0]
        # get proteins and sort by descending length
        gene_prots = v
        sorted_prots = sorted(list(gene_prots.keys()), key=len, reverse=True)

        protid = 1
        alleleid = 1
        total_alleles = 0
        dna_file = '{0}_{1}'.format('dna', g)
        dna_file = os.path.join(output_dir, dna_file)
        protein_file = '{0}_{1}'.format('protein', g)
        protein_file = os.path.join(output_dir, protein_file)
        for prot in sorted_prots:
            current_prot = gene_prots[prot]

            alleles = [a[0] for a in current_prot[0]]
            total_alleles += len(alleles)
            sequences = [a[1] for a in current_prot[0]]

            prot_header = '>{0}*{1}\n'.format(gene_id, protid)
            with open(protein_file, 'a') as pfile:
                prot_sequence = prot
                prot_record = '{0}{1}\n'.format(prot_header,
                                                prot_sequence)
                pfile.write(prot_record)

            with open(dna_file, 'a') as dfile:
                for i, allele in enumerate(alleles):
                    dna_header = '>{0}*{1}.{2}_{3}\n'.format(gene_id,
                                                             protid,
                                                             alleleid,
                                                             allele)
                    dna_sequence = sequences[i]
                    dna_record = '{0}{1}\n'.format(dna_header,
                                                   dna_sequence)
                    dfile.write(dna_record)
                    alleleid += 1

            alleleid = 1
            protid += 1

    return [gene_id, len(sorted_prots), total_alleles]


def protein_diversity(file, output_dir):
    """ Determines DNA sequences that code for the same protein
        and writes new files with the DNA sequences and protein
        sequences, including dot notation in the DNA sequences
        headers to highlight the protein they code and if there
        are more DNA sequences coding for the same  protein.

        Args:
            file (str): path a FASTA file with DNA sequences
            for alleles of a single gene.
            output_dir (str): output directory that will be used
            to store the new files.

        Returns:
            Writes a FASTA file with the DNA sequences for
            a gene, including dot notation in the header and
            a FASTA file with the distinct protein sequences
            resulting from the translation of all DNA sequences,
            ordered by descending length.
    """

    try:
        gene_proteins = group_by_protein(file)
        res = write_files(gene_proteins, output_dir)
        print('Processed file: {0}'.format(file))
        return res
    except Exception:
        print('Failed for file: {0}'.format(file))


def plot_diversity(data, title, output_dir):
    """ Creates a HTML file that can be opened with a Browser
        to visualize a Bar plot with the number of distinct
        proteins and distinct alleles per schema gene.

        Args:
            data (list): list with sublists. Each sublist has
            3 variables: gene identifier, number of distinct
            proteins and number of distinct alleles.
            title (str): title for the Bar plot.
            output_dir (str): output diretory where the HTML
            file will be created.

        Returns:
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
    fig.update_layout(barmode='group', title_text=title)

    output_file = os.path.join(output_dir, '{0}.html'.format(title))

    plot(fig, filename=output_file, auto_open=False)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True, dest='input_files',
                        help='The chewBBACA schema directory.')

    parser.add_argument('-o', type=str, required=True, dest='output_directory',
                        help='Output directory that will be created to store'
                        ' the new files.')

    parser.add_argument('-t', type=int, required=False, dest='cpu_cores',
                        default=1, help='Number of CPU cores to use.')

    parser.add_argument('-p', type=str, required=False, dest='create_plot',
                        choices=['no', 'yes'], default='no', help='If the'
                        ' process should create a Bar plot comparing the number'
                        ' of alleles and the number of distinct proteins for each'
                        ' gene.')

    args = parser.parse_args()

    return [args.input_files, args.output_directory, args.cpu_cores, args.create_plot]


def main(schema_dir, output_dir, cpu_cores, create_plot):

    # get list of schema files
    schema_files = [os.path.join(schema_dir, file) for file in os.listdir(schema_dir) if '.fasta' in file]
    schema_basename = os.path.basename(schema_dir) if schema_dir[-1] != '/' else os.path.basename(schema_dir[:-1])

    # create directory to store new files
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # cluster alleles identifiers based on protein sequence
    plot_data = []
    pool = multiprocessing.Pool(cpu_cores)
    for file in schema_files:
        pool.apply_async(protein_diversity, (file, output_dir), callback=plot_data.append)
    pool.close()
    pool.join()

    if create_plot == 'yes':
        plot_diversity(plot_data, schema_basename, output_dir)

    print('Done!')


if __name__ == '__main__':

    args = parse_arguments()
    main(args[0], args[1], args[2], args[3])
