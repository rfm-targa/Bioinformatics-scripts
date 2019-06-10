#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR
    
    Rafael Mamede

DESCRIPTION
        
    Creates a basic Heatmap from a Snippy (created by Torsten Seemann and available at 
    https://github.com/tseemann/snippy) core.aln or core.full.aln file.
    
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
import itertools

import plotly.offline
from Bio import SeqIO
from Bio.Alphabet import generic_dna


#snippy_aln_file = '/home/rfm/Desktop/rfm/Lab_Analyses/GBS_all_samples_reads/GBSAC13832_1_ref/core.full.aln'
#reference_id = 'gbs'
#plot_filename = 'gbs'


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
    snp_hoverinfo = []
    
    # dictionary to change each sequence character to an integer
    differences_mapping = {'A':1, 'T':1, 'G':1, 'C':1,
                           'a':1, 't':1, 'g':1, 'c':1, 
                           '-':2, 'N':3, 'X':3, 'n':3}
    
    color_range = []
    for assembly, snp_profile in assemblies_profiles.items():
        
        # sample identifier for yaxis labels
        yaxis_labels.append(assembly)
        
        # sample integer vectors, new list has an integer for each character in the original sequence
        reference_based = [0 if snp_profile[i] == assemblies_profiles[reference][i] \
                           else differences_mapping[snp_profile[i]] for i in range(len(assemblies_profiles[reference]))]
        snp_tracers_data.append(reference_based)
        
        # append distinct integers to know which type of chars are present in the alignment
        color_range.append(list(set(reference_based)))
        
        # create hovertext to display nucleotides or special characters while hovering over plot
        sample_hovertext = [snp_profile[i] if snp_profile[i] == assemblies_profiles[reference][i] \
                            else '{0} ({1})'.format(snp_profile[i], assemblies_profiles[reference][i]) for i in range(len(snp_profile))]
        snp_hoverinfo.append(sample_hovertext)


    print('Generating Heatmap tracer...')
    
    # colors for each type of char
    colors = {0: '#ece7f2', 1: '#41ab5d',
              2: '#0868ac', 3: '#252525'}
    
    # flatten and keep distinct group of integers
    color_range = list(set(list(itertools.chain.from_iterable(color_range))))
    # sort integers list
    color_range.sort()
    
    # find max integer and determine the number of levels for the colorscale
    max_val = max(color_range)
    normalized_color_range = [c/max_val for c in color_range]
    
    # create colorscale
    colorscale = []
    for i in range(len(color_range)):
        colorscale.append([normalized_color_range[i], colors[color_range[i]]])

    # create Heatmap tracer
    # created as normal Python dictionary for speed
    # go.Heatmap takes a long time with large arrays
    snp_heatmap_trace = dict(z = snp_tracers_data,
                             colorscale = colorscale,
                             y = yaxis_labels,
                             showscale = False,
                             type = 'heatmap',
                             # creating hovertext with 'x', 'y' and 'z' info can
                             # lead to a big html file
                             hoverinfo= 'x+y+text',
                             text = snp_hoverinfo,
                             # add gap between each sample row
                             ygap =	3)
    
    snp_heatmap_data = [snp_heatmap_trace]
    
    # integer to feature mapping
    feature_mapping = {0: 'Equal', 1: 'SNP',
                       2: 'Zero coverage or deletion (-)', 
                       3: 'Low coverage (N)<br>Masked region on reference (X)<br>Heterozygous or poor quality genotype (n)'}
    
    # create ghost tracer to have custom legend
    # create a tracer without data but with a legend to show the meanign of each
    # heatmap color
    for f in range(len(color_range)):
        snp_ghost_trace = dict(x = [None],
                               y = [None],
                               mode = 'markers',
                               marker = dict(size=10, color=colors[color_range[f]], symbol = 'square',
                                             line = dict(color = '#000000', width = 1)),
                               showlegend = True,
                               name = feature_mapping[color_range[f]])
    
        snp_heatmap_data.append(snp_ghost_trace)
    
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
    plotly.offline.plot(snp_heatmap_fig, filename = output_plot_file, auto_open = False, validate=False)

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

