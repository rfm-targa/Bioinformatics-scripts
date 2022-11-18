#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
Creates a basic Heatmap from a Snippy (https://github.com/tseemann/snippy)
core.aln or core.full.aln file.

Each sample nucleotide position is colored according to the same position
in the reference sequence.
Colorscale hex codes:
Equal to reference: #ece7f2
SNPs: #41ab5d
Zero coverage or deletions: #0868ac
Low coverage, masked region on reference and heterozygous or poor quality
genotype: #252525

This script requires that 'Biopython' and 'Plotly' be installed within the
Python 3 environment that is used to run the script.

Note: the output file size will vary according to the number of samples.
When using a core.full.aln file, the output HTML file might be large
depending on the number of samples (e.g: a core.full.aln file with the
alignment between the assembly of 1 reference and 7 samples of Streptococcus
agalactiae originated a HTML file with 124MB.
"""


import argparse
import itertools

import plotly.offline
from Bio import SeqIO


def heatmap_trace(y_data, z_data, text, colorscale='Viridis',
                  showscale=False, hoverinfo='x+y+text', sample_gap=3):
    """Create a dictionary with data to create a Heatmap trace.

    Parameters
    ----------
    y_data : list
        Sample labels for the yaxis.
    z_data : list
        Values for the zaxis, one sublist of values for each sample.
    text : list
        Text info that will be displayed for the zaxis hoverinfo.
    colorscale : list
        Sublists with one integer value for the scale value and a
        string with the hexcode for the color (default Viridis).
    showscale : bool
        Display heatmap vertical colorscale.
    hoverinfo : str
        The info that is displayed while hovering through the plot.
    sample_gap : int
        The space between each sample row in the heatmap.

    Returns
    -------
    heatmap_trace : dict
        Dictionary with keys and values for the heatmap trace.
    """
    heatmap_trace = dict(z=z_data,
                         colorscale=colorscale,
                         y=y_data,
                         showscale=showscale,
                         # creating custom hovertext for 'x', 'y' and 'z'
                         # info can lead to a big html file
                         hoverinfo=hoverinfo,
                         text=text,
                         # add gap between each sample row
                         ygap=sample_gap)

    return heatmap_trace


def legend_ghost_tracer(trace_name, x_data=[None], y_data=[None],
                        trace_mode='markers', marker_size=10,
                        marker_color='#000000', marker_symbol='square',
                        marker_line_color='#000000', marker_line_width=1):
    """Create a dictionary with data to create an invisible legend trace.

    Parameters
    ----------
    trace_name : str
        Name displayed in the legend.
    x_data : list
        xaxis data to display in plot.
    y_data : list
        yaxis data to display in plot.
    trace_mode : str
        Display mode of data points.
    marker_size : int
        Size of the markers.
    marker_color : str
        Color filling the marker symbol.
    marker_symbol : str
        Shape of the marker symbol.
    marker_line_color : str
        Color of the line highlighting the marker symbol.
    marker_line_width : int
        Width of the line highlighting the marker symbol.

    Returns
    -------
    ghost_trace : dict
        Dictionary with the data to create an invisible scatter
        that displays hover info.
    """
    ghost_tracer = dict(x=x_data,
                        y=y_data,
                        mode=trace_mode,
                        marker=dict(size=marker_size, color=marker_color,
                                    symbol=marker_symbol,
                                    line=dict(color=marker_line_color,
                                              width=marker_line_width)),
                        showlegend=True,
                        name=trace_name)

    return ghost_tracer


def figure_layout(plot_title, xaxis_title, yaxis_title):
    """Create the dictionary with data for the Figure object layout.

    Parameters
    ----------
    plot_title : str
        Main title for the plot.
    xaxis_title : str
        Title of the xaxis.
    yaxis_title : str
        Title of the yaxis.

    Returns
    -------
    layout : dict
        Dictionary with keys and values for the layout features.
    """
    layout = dict(title=plot_title,
                  xaxis=dict(title=xaxis_title),
                  yaxis=dict(title=yaxis_title),
                  autosize=True,
                  margin=dict(l=150, r=120, b=120)
                  )

    return layout


def main(input_file, output_file, reference_id):

    # passed arguments
    snp_align_file = input_file
    # custom identifier for reference
    reference = '{0} (ref)'.format(reference_id)
    output_plot_file = '{0}.html'.format(output_file)

    print('\nImporting and processing alignment data...')

    # read SNPs profiles into dictionary
    reference_count = 0
    sample_count = 0
    assemblies_profiles = {}
    for snp_profile in SeqIO.parse(snp_align_file, 'fasta'):
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
        print('Found {0} reference and {1} samples.'.format(reference_count,
                                                            sample_count))
    else:
        raise Exception('No reference found in the Snippy alignment file! '
                        'Please make sure that the .aln file has a sequence '
                        'with ">Reference" as header identifier.')

    print('Determining differences between reference and remaining samples...')

    # process SNP data to get tracers for heatmap
    yaxis_labels = []
    snp_tracers_data = []
    snp_hoverinfo = []

    # dictionary to change each sequence character to an integer
    differences_mapping = {'A': 1, 'T': 1, 'G': 1, 'C': 1,
                           'a': 1, 't': 1, 'g': 1, 'c': 1,
                           '-': 2, 'N': 3, 'X': 3, 'n': 3}

    color_range = []
    for assembly, snp_profile in assemblies_profiles.items():
        # sample identifier for yaxis labels
        yaxis_labels.append(assembly)

        # sample integer vectors, new list has an integer for each
        # character in the original sequence
        reference_based = [0 if snp_profile[i] == assemblies_profiles[reference][i] \
                           else differences_mapping[snp_profile[i]]
                           for i in range(len(assemblies_profiles[reference]))]
        snp_tracers_data.append(reference_based)

        # append distinct integers to know which type of chars are
        # present in the alignment
        color_range.append(list(set(reference_based)))

        # create hovertext to display nucleotides or special characters
        # while hovering over plot
        sample_hovertext = [snp_profile[i]
                            if snp_profile[i] == assemblies_profiles[reference][i] \
                            else '{0} ({1})'.format(snp_profile[i], assemblies_profiles[reference][i])
                            for i in range(len(snp_profile))]
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

    # create Heatmap trace
    # created as normal Python dictionary for speed
    # go.Heatmap takes a long time with large arrays
    snp_heatmap_trace = heatmap_trace(yaxis_labels, snp_tracers_data,
                                      snp_hoverinfo, colorscale=colorscale)

    snp_heatmap_data = [snp_heatmap_trace]

    # integer to feature mapping
    feature_mapping = {0: 'Equal', 1: 'SNP',
                       2: 'Zero coverage or deletion (-)',
                       3: 'Low coverage (N)<br>Masked region on reference '
                          '(X)<br>Heterozygous or poor quality genotype (n)'}

    # create ghost trace for custom legend
    # create a trace without data but with a legend to show the meaning of each
    # heatmap color
    for f in range(len(color_range)):
        legend_tracer = legend_ghost_tracer(feature_mapping[color_range[f]],
                                            marker_color=colors[color_range[f]])

        snp_heatmap_data.append(legend_tracer)

    # figure layout
    snp_heatmap_layout = figure_layout('{0} Heatmap'.format(reference),
                                       'Nucleotide Position',
                                       'Sample Identifier')

    print('Creating Figure object and Heatmap HTML file...')

    # create Figure, also as a Python dictionary for speed
    snp_heatmap_fig = dict(data=snp_heatmap_data, layout=snp_heatmap_layout)
    # plot the Heatmap, no validation to be faster
    plotly.offline.plot(snp_heatmap_fig, filename=output_plot_file,
                        auto_open=False, validate=False)

    print('Done!\n')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str, required=True,
                        dest='input_file',
                        help='Path to the Snippy multi-fasta alignment file. '
                             'Either the core.aln with concatenated SNPs or '
                             'the core.full.aln with the full alignment.')

    parser.add_argument('-o', '--output-file', type=str, required=True,
                        dest='output_file',
                        help='Path to the output HTML file with the heatmap '
                             '(.html extension will be added to end of the '
                             'filename).')

    parser.add_argument('-rid', '--reference-id', type=str, required=True,
                        dest='reference_id',
                        help='User-defined identifier string for the '
                             'reference sample.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
