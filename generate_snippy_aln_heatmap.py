#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR
    
    Rafael Mamede

DESCRIPTION
        
    Creates a basic Heatmap from a Snippy core.aln or core.full.aln file.
    
    Each nucleotide site is colored according to the reference sequence.
    Sites with the same nucleotide as the reference are colored in gray and
    differences are colored in blue.
    
    This script requires that 'Biopython' and 'Plotly' be installed within the 
    Python 3 environment that is used to run this script.

    This file can also be imported as a module and contains the following
    functions:

        * snippy_heatmap - creates a Heatmap to highlight nucleotide differences
          between a reference and other samples;
        * main - the main function of the script, used to parse arguments from 
          the command line and call the snippy_heatmap function.
"""


import argparse

import plotly.offline
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def snippy_heatmap(snippy_aln_file, reference_id, plot_filename):
    """ Creates a Plotly Heatmap based on sequence differences present in the 
        Snippy .aln files, either the core.aln with concatenated SNPs or the 
        core.full.aln file with full assembly/genome multi-alignment.
    
        Args: 
            snippy_aln_file (str): path to the Snippy .aln file.
            reference_id (str): identifier for the 'Reference' sequence in the
            Snippy .aln file.
            plot_filename (str): prefix for the HTML file that will be used to
            visualize the Heatmap.
        
        Returns:
            Creates a 'plot_filename.html' file in the working directory.
    """

    # passed arguments
    snp_align_file = snippy_aln_file
    reference = '{0} (ref)'.format(reference_id)
    output_plot_file = '{0}.html'.format(plot_filename)

    print('\nImporting and processing alignment data...')
    
    # read SNPs profiles into dictionary
    reference_count = 0
    sample_count = 0
    assemblies_profiles = {}
    for snp_profile in SeqIO.parse(snp_align_file, 'fasta', generic_dna):
        assembly_id = snp_profile.id
        assembly_profile = str(snp_profile.seq)

        if 'Reference' in assembly_id:
            assemblies_profiles[reference] = assembly_profile
            reference_count += 1
        else:
            assemblies_profiles[assembly_id] = assembly_profile
            sample_count += 1
    
    if reference_count > 0:
        print('Found {0} reference and {1} samples.'.format(reference_count, sample_count))
    else:
        raise Exception('No reference found in the Snippy alignment file! Please make '
                        'sure that the .aln file has a sequence with ">Reference" as header identifier.')

    print('Determining differences between reference and remaining samples...')

    # process SNP data to get tracers for heatmap
    yaxis_labels = []
    snp_tracers_data = []

    for assembly, snp_profile in assemblies_profiles.items():
        
        yaxis_labels.append(assembly)
        reference_based = ['0' if snp_profile[i] == assemblies_profiles[reference][i] \
                           else '1' for i in range(len(assemblies_profiles[reference]))]
        
        snp_tracers_data.append(reference_based)

    print('Generating Heatmap tracer...')

    # create Heatmap tracer
    # created as normal Python dictionary for speed
    # go.Heatmap takes a long time with large arrays
    snp_heatmap_trace = dict(y = yaxis_labels,
                             z = snp_tracers_data,
                             colorscale = [[0, 'rgb(224, 224, 224)'],
                                           [1, 'rgb(51,153,255)']],
                             showscale = False,
                             type = 'heatmap')
    
    snp_heatmap_data = [snp_heatmap_trace]
    
    snp_heatmap_layout = dict(title = '{0} Heatmap'.format(reference),
                              xaxis = dict(title = 'Nucleotide Position'),
                              yaxis = dict(title = 'Sample Identifier'),
                              autosize = True,
                              margin = dict(l = 150, r = 120, b = 120)
                              )
    
    print('Creating Figure object and Heatmap HTML file...')
    
    # create Fig, also as a Python dictionary for speed
    snp_heatmap_fig = dict(data = snp_heatmap_data, layout = snp_heatmap_layout)
    # plot the Heatmap, no validation to be faster
    plotly.offline.plot(snp_heatmap_fig, filename = output_plot_file, auto_open = False, 
                        validate=False)

    print('Done!\n')


def main():
    
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument('snippy_aln_file', type=str, 
                        help='Path to the Snippy multi-fasta alignment file.'
                             ' Either the core.aln with concatenated SNPs or'
                             ' the core.full.aln with full alignment.')
    
    parser.add_argument('reference_id', type=str,
                       help='User-defined identifier string for the reference sample.')
    
    parser.add_argument('plot_filename', type=str,
                        help='Name for the output HTML file with the heatmap'
                             ' (.html extension will be added to end of filename.)')
    
    args = parser.parse_args()
    snippy_heatmap(args.snippy_aln_file, args.reference_id, args.plot_filename)
    

if __name__ == '__main__':
    
    main()

