#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR
    
    Rafael Mamede
    github: @rfm-targa

DESCRIPTION
        
    Creates a basic Heatmap from a Snippy (created by Torsten Seemann and available at 
    https://github.com/tseemann/snippy) core.aln or core.full.aln file.
    
    Each sample nucleotide position is colored according to the same position in the 
    reference sequence.
    Colorscale hex codes:
    Equal to reference: #ece7f2
    SNPs: #41ab5d
    Zero coverage or deletions: #0868ac
    Low coverage, masked region on reference and heterozygous or poor quality genotype: #252525
    
    This script requires that 'Biopython' and 'Plotly' be installed within the 
    Python 3 environment that is used to run this script.

    This file can also be imported as a module and contains the following
    functions:
        
        * create_heatmap_tracer - creates a tracer for a heatmap plot;
        * legend_ghost_tracer - creates a tracer with no data but that displays a
                                marker for a custom legend;
        * figure_layout - creates the layout for the plot figure;
        * snippy_heatmap - creates a Heatmap plot from a Snippy .aln file to 
                           highlight sequence differences between a reference 
                           and other samples;
        * main - the main function of the script, used to parse arguments from 
          the command line and call the snippy_heatmap function.
    
    Note: the output file size will vary according to the number of samples. When using 
    a core.full.aln file, the output HTML file might be large depending on the number of samples.
    e.g: a core.full.aln file with the alignment between the assembly of 1 reference and 7 samples
    of Streptococcus agalactiae originated a HTML file with 124MB.
"""


import argparse
import itertools

import plotly.offline
from Bio import SeqIO
from Bio.Alphabet import generic_dna


def create_heatmap_tracer(y_data, z_data, text, colorscale='Viridis', showscale=False, 
                          plot_type='heatmap', hoverinfo='x+y+text', sample_gap=3):
    """ Creates a tracer for a heatmap plot.
        
        Args: 
            y_data (list): sample labels for the yaxis;
            z_data (list): values for the zaxis, one sublist of values for each 
                           sample;
            text (list): text info that will be displayed for the zaxis hoverinfo;
            colorscale (list): sublists with one integer value for the scale value
                               and a string with the hexcode for the color (default Viridis);
            showscale (bool): display heatmap vertical colorscale (default False);
            plot_type (str): type of plot (default heatmap);
            hoverinfo (str): the info that is displayed while hovering through the 
                             plot (default x+y+text);
            sample_gap (int): the space between each sample row in the heatmap (default 3).
        
        Returns:
            heatmap_tracer (dict): keys and values for each data feature of the trace.
    """

    heatmap_tracer = dict(z = z_data,
                         colorscale = colorscale,
                         y = y_data,
                         showscale = showscale,
                         type = plot_type,
                         # creating custom hovertext for 'x', 'y' and 'z' info can
                         # lead to a big html file
                         hoverinfo = hoverinfo,
                         text = text,
                         # add gap between each sample row
                         ygap = sample_gap)

    return heatmap_tracer


def legend_ghost_tracer(tracer_name, x_data=[None], y_data=[None], tracer_mode='markers', 
                        marker_size=10, marker_color='#000000', marker_symbol='square', 
                        marker_line_color='#000000', marker_line_width=1):
    """ Creates a tracer for a heatmap plot.
        
        Args: 
            tracer_name (str): name displayed in the legend;
            x_data (list): xaxis data to display in plot (default None);
            y_data (list): yaxis data to display in plot (default None);
            tracer_mode (str): display mode of data points (default markers);
            marker_size (int): size of the markers (default 10);
            marker_color (str): color filling the marker symbol (default black);
            marker_symbol (str): shape of the marker symbol (default square);
            marker_line_color (str): color of the line highlighting the marker symbol (default black);
            marker_line_width (int): width of the line highlighting the marker symbol (default black).
        
        Returns:
            ghost_tracer (dict): keys and values for each data feature of the trace.
    """

    ghost_tracer = dict(x = x_data,
                        y = y_data,
                        mode = tracer_mode,
                        marker = dict(size=marker_size, color=marker_color, 
                                      symbol=marker_symbol,
                                      line = dict(color=marker_line_color, 
                                                  width=marker_line_width)),
                        showlegend = True,
                        name = tracer_name)

    return ghost_tracer


def figure_layout(plot_title, xaxis_title, yaxis_title):
    """ Creates a tracer for a heatmap plot.
        
        Args: 
            plot_title (str): main title of the plot;
            xaxis_title (str): title of the xaxis;
            yaxis_title (str): title of the yaxis.
        
        Returns:
            layout (dict): dictionary with keys and values for the layout features.
    """

    layout = dict(title = plot_title,
                  xaxis = dict(title=xaxis_title),
                  yaxis = dict(title=yaxis_title),
                  autosize = True,
                  margin = dict(l=150, r=120, b=120)
                  )
    
    return layout


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
    # custom identifier for reference
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
    
    # proceed only if there is a reference
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
    snp_heatmap_trace = create_heatmap_tracer(yaxis_labels, snp_tracers_data, snp_hoverinfo, colorscale=colorscale)

    snp_heatmap_data = [snp_heatmap_trace]

    # integer to feature mapping
    feature_mapping = {0: 'Equal', 1: 'SNP',
                       2: 'Zero coverage or deletion (-)', 
                       3: 'Low coverage (N)<br>Masked region on reference (X)<br>Heterozygous or poor quality genotype (n)'}

    # create ghost tracer for custom legend
    # create a tracer without data but with a legend to show the meaning of each
    # heatmap color
    for f in range(len(color_range)):
        legend_tracer = legend_ghost_tracer(feature_mapping[color_range[f]], marker_color = colors[color_range[f]])
    
        snp_heatmap_data.append(legend_tracer)

    # figure layout
    snp_heatmap_layout = figure_layout('{0} Heatmap'.format(reference), 'Nucleotide Position', 'Sample Identifier')

    print('Creating Figure object and Heatmap HTML file...')

    # create Figure, also as a Python dictionary for speed
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
                             ' the core.full.aln with the full alignment.')
    
    parser.add_argument('reference_id', type=str,
                       help='User-defined identifier string for the reference sample.')
    
    parser.add_argument('plot_filename', type=str,
                        help='Name for the output HTML file with the heatmap'
                             ' (.html extension will be added to end of filename).')
    
    args = parser.parse_args()
    snippy_heatmap(args.snippy_aln_file, args.reference_id, args.plot_filename)
    

if __name__ == '__main__':
    
    main()

