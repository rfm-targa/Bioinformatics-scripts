#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script computes intracluster and intercluster distances
based on predefined clusters and allele calling results. The
script creates several files per cluster with statistics for
the wgMLST and the cgMLST and HTML files with interactive plots
(boxplots with the distribution of distance values, heatmaps
representing the distance matrix and dendograms representing
the results of hierarchical clustering).

Code documentation
------------------
"""


import os
import io
import csv
import shutil
import argparse
import contextlib

import numpy as np
import pandas as pd
from plotly.offline import plot
import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.cluster.hierarchy import dendrogram, linkage

# these modules must be in the same directory
import mask_matrix as mm
import distance_matrix as dm
import Extract_cgAlleles as cg


def nostdout(function, args):
    """Silence stdout from a function."""
    with contextlib.redirect_stdout(io.StringIO()):
        return function(*args)


def simply_return(data):
    """Return provided argument."""
    return data


def read_tabular(input_file, delimiter='\t'):
    """Read a tabular file.

    Parameters
    ----------
    input_file : str
        Path to a tabular file.
    delimiter : str
        Delimiter used to separate file fields.

    Returns
    -------
    lines : list
        A list with a sublist per line in the input file.
        Each sublist has the fields that were separated by
        the defined delimiter.
    """
    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines


def intracluster_stats(distance_matrix, row_id):
    """Compute the min, max and mean distances to samples in same cluster.

    Parameters
    ----------
    distance_matrix : pandas.core.frame.DataFrame
        Pandas dataframe that represents a distance matrix
        with the number of allelic differences.
    row_id : str
        Sample identifier.

    Returns
    -------
    sample_stats : list
        A list with the minimum, maximum and mean distance
        between the sample with provided identifier and
        samples in the same cluster.
    """
    # get row
    current_row = distance_matrix.loc[row_id]
    # drop self
    current_row = current_row.drop([row_id])
    # get minimum distance in same cluster
    minimum_distance = current_row.nsmallest(1, keep='all')
    stats_min = (list(minimum_distance.index),
                 minimum_distance[0])
    # get maximum distance in same cluster
    maximum_distance = current_row.nlargest(1, keep='all')
    stats_max = (list(maximum_distance.index),
                 maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = current_row.mean()
    stats_mean = mean_distance

    sample_stats = [stats_min, stats_max, stats_mean]

    return sample_stats


def intercluster_stats(distance_matrix, row_id):
    """Compute the min, max and mean distances to samples in other clusters.

    Parameters
    ----------
    distance_matrix : pandas.core.frame.DataFrame
        Pandas dataframe that represents a distance matrix
        with the number of allelic differences.
    row_id : str
        Sample identifier.

    Returns
    -------
    sample_stats : list
        A list with the minimum, maximum and mean distance
        between the sample with provided identifier and
        samples in other clusters.
    """
    # get minimum distance to strains in other clusters
    minimum_distance = distance_matrix[row_id].nsmallest(1, keep='all')
    stats_min = (list(minimum_distance.index),
                 minimum_distance[0])
    # get maximum distance to strains in other clusters
    maximum_distance = distance_matrix[row_id].nlargest(1, keep='all')
    stats_max = (list(maximum_distance.index),
                 maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = distance_matrix[row_id].mean()
    stats_mean = mean_distance

    sample_stats = [stats_min, stats_max, stats_mean]

    return sample_stats


def select_centroids(cluster_stats):
    """Select cluster centroids based on the mean distance between all samples.

    Parameters
    ----------
    cluster_stats : dict
        Dictionary with distance statistics for all
        samples in a cluster.

    Returns
    -------
    centroid : str
        Identifier of the sample selected as centroid.
    """
    if len(cluster_stats) > 1:
        means = [(i, j['mean']) for i, j in cluster_stats.items()]
        # sort based on decreasing mean distance
        sorted_means = sorted(means, key=lambda x: x[1])
        centroid = sorted_means[0][0]
    # select singleton as centroid
    else:
        centroid = list(cluster_stats.keys())[0]

    return centroid


def centroids_inter_dists(distance_matrix, centroids_ids, identifier):
    """Compute distances between cluster centroids.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.
    centroids_ids : list
        List with the identifiers for all centroids.
    identifier : str
        Sample identifier for the current centroid.

    Returns
    -------
    centroids_stats : list
        List with the minimum, maximum and mean distance
        to centroids from other clusters.
    """
    current_cols = ['FILE', identifier]
    df = pd.read_csv(distance_matrix, sep='\t',
                     usecols=current_cols, index_col='FILE')
    # exclude self and all non-centroids
    outdf = df[df.index.isin(centroids_ids)]
    # get minimum distance to other centroids
    minimum_distance = outdf[identifier].nsmallest(1, keep='all')
    centroid_min = (list(minimum_distance.index),
                    minimum_distance[0])
    # get maximum distance to strains in other clusters
    maximum_distance = outdf[identifier].nlargest(1, keep='all')
    centroid_max = (list(maximum_distance.index),
                    maximum_distance[0])
    # get mean distance in same cluster
    mean_distance = outdf[identifier].mean()
    centroid_mean = mean_distance

    centroids_stats = [centroid_min, centroid_max, centroid_mean]

    return centroids_stats


def compute_intracluster_stats(distance_matrix, cluster_ids):
    """Compute the min, max and mean distance between all samples in a cluster.

    Also selects the cluster centroid based on the minimum mean
    intracluster distance.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.
    cluster_ids : list
        List with the identifiers of the samples in the
        cluster.

    Returns
    -------
    cluster_stats : dict
        Dictionary with sample identifiers as keys and
        a dictionary with key:value pairs for the minimum,
        maximum and mean distance values. The dictionary
        also includes a key:value pair for the selected
        centroid.
    """
    # only read columns for the samples that belong to current cluster
    current_cols = ['FILE'] + cluster_ids
    df = pd.read_csv(distance_matrix, sep='\t',
                     usecols=current_cols, index_col='FILE')

    cluster_stats = {e: {} for e in cluster_ids}

    # do not compute statistics for singletons
    if len(cluster_ids) > 1:
        # compute statistics for the cluster
        for e in cluster_ids:
            strain_stats = intracluster_stats(df, e)
            cluster_stats[e]['min'] = strain_stats[0]
            cluster_stats[e]['max'] = strain_stats[1]
            cluster_stats[e]['mean'] = strain_stats[2]
    # add empty strings to singletons stats
    else:
        for e in cluster_ids:
            cluster_stats[e]['min'] = ([''], '')
            cluster_stats[e]['max'] = ([''], '')
            cluster_stats[e]['mean'] = ''

    # determine centroid candidates based on minimum mean
    # intracluster distance
    cluster_centroid = select_centroids(cluster_stats)
    cluster_stats['centroids'] = cluster_centroid

    return cluster_stats


def clusters_stats(distance_matrix, clusters):
    """Compute intracluser stats for all clusters.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.
    clusters: dict
        Dictionary with cluster identifiers as keys and
        a list with the identifiers of samples in each
        cluster as values.

    Returns
    -------
    stats : dict
        Dictionary with cluster identifiers as keys and
        a dictionary with the intracluster stats for all
        samples in each cluster as values.
    """
    # use Pandas to get columns for each cluster
    stats = {}
    for k, v in clusters.items():
        current_stats = compute_intracluster_stats(distance_matrix, v)
        stats[k] = current_stats

    return stats


def clusters_boxplot(distance_matrix, clusters):
    """Create Boxplots to plot the distribution of intracluster distances.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.
    clusters : dict
        Dictionary with cluster identifiers as keys and
        a list with the identifiers of samples in each
        cluster as values.

    Returns
    -------
    traces : list
        List with one Boxplot trace per cluster.
    """
    # use Pandas to get columns for each cluster
    traces = []
    for k, v in clusters.items():
        current_cols = ['FILE'] + v
        df = pd.read_csv(distance_matrix, sep='\t',
                         usecols=current_cols, index_col='FILE')

        # get lines with cluster samples
        df = df.loc[current_cols[1:]]

        # get distance values
        cluster_array = df.to_numpy()
        upper_values = cluster_array[np.triu_indices(len(v), k=1)]
        values = upper_values.tolist()

        if len(values) >= 5:
            # create trace for cluster's boxplot
            trace = go.Box(y=values, name=k,
                           marker=dict(color='#2b8cbe',
                                       line=dict(width=1, color='#252525')),
                           fillcolor='#2b8cbe',
                           line_color='#252525',
                           boxpoints='outliers', jitter=0.5)
        else:
            # only display points if n<5
            trace = go.Box(y=values, name=k,
                           marker=dict(color='#2b8cbe',
                                       line=dict(width=1, color='#252525')),
                           fillcolor='rgba(0,0,0,0)',
                           line_color='rgba(0,0,0,0)',
                           pointpos = 0,
                           boxpoints='all',
                           jitter=0.5)

        traces.append(trace)

    return traces


def clusters_heatmap(distance_matrix):
    """Create a Heatmap represeting a distance matrix.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.

    Returns
    -------
    heatmap_trace : go.Heatmap
        Heatmap trace to represent the distance matrix
        of allelic distances for the cluster.
    """
    df = pd.read_csv(distance_matrix,
                     sep='\t',
                     index_col='FILE')

    # get distance values
    cluster_array = df.to_numpy()
    values = cluster_array.tolist()

    matrix_header = list(df.columns)

    # invert nested list with values to represent
    # the distance matrix as it is in the TSV file
    heatmap_trace = go.Heatmap(z=values[::-1],
                               x=matrix_header,
                               y=matrix_header[::-1],
                               colorscale='Viridis')

    return heatmap_trace


def cluster_dendogram(distance_matrix, output_file,
                      linkage_function, distance_function):
    """Clusters samples based on a linkage function and a distance function.

    Parameters
    ----------
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.
    output_file : str
        Path to the output file.
    linkage_function : func
        Function to use to perform hierarchical clustering.
    distance_function : func
        Function to use to compute the distance matrix.

    Returns
    -------
    True
    """
    # Read distance matrix
    distance_df = pd.read_csv(distance_matrix,
                              sep='\t',
                              index_col='FILE')

    labels = list(distance_df.index)
    # do not create dendogram if we only have one sample
    if len(labels) < 2:
        return False

    # Dataframe to Numpy array
    distance_array = distance_df.to_numpy()
    # ff.create_dendogram only accepts condensed matrices
    condensed_array = np.triu(distance_array)

    # Initialize figure by creating upper dendrogram
    fig = ff.create_dendrogram(condensed_array, orientation='bottom',
                               labels=labels, distfun=distance_function,
                               linkagefun=linkage_function)

    # hide upper dendogram
    fig.for_each_trace(lambda trace: trace.update(visible=False))

    # get ordered sample labels from upper dendogram
    ordered_labels = list(fig.layout['xaxis']['ticktext'])

    # change yaxis for dendogram traces
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(condensed_array, orientation='right',
                                       distfun=distance_function,
                                       linkagefun=linkage_function)

    # change xaxis for side dendogram
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data)

    # Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(map(int, dendro_leaves))
    heat_data = distance_array
    # reorder data points so that they align with the dendogram leaves
    # align with side dendogram
    heat_data = heat_data[dendro_leaves,:]
    # align with upper dendogram
    heat_data = heat_data[:,dendro_leaves]

    # create list with hovertext for each heatmap cell
    heatmap_hovertext = []
    for i, l in enumerate(ordered_labels):
        heatmap_hovertext.append([])
        for j, l2 in enumerate(ordered_labels):
            heatmap_hovertext[-1].append('x: {0}<br>y: {1}<br>Distance: {2}'.format(l, l2, heat_data[i][j]))

    heatmap = [go.Heatmap(x=dendro_leaves,
                          y=dendro_leaves,
                          z=heat_data,
                          colorscale='Blues',
                          hoverinfo='text',
                          text=heatmap_hovertext)]

    # adjust heatmap tick values position
    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data)

    # Edit Layout
    fig.update_layout({'width': 1200, 'height': 900,
                       'showlegend': False, 'hovermode': 'closest'})

    # Edit xaxis with heatmap
    fig.update_layout(xaxis={'domain': [.15, 1],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'showticklabels': False,
                             'ticks': ''})

    # Edit xaxis2 with dendogram
    fig.update_layout(xaxis2={'domain': [0, .15],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ''})

    # Edit yaxis with heatmap and side dendogram
    fig.update_layout(yaxis={'domain': [0, 1],
                             'mirror': False,
                             'showgrid': False,
                             'showline': False,
                             'zeroline': False,
                             'showticklabels': False,
                             'ticks': ''})

    # Edit yaxis2 with top dendogram
    fig.update_layout(yaxis2={'domain': [.825, .975],
                              'mirror': False,
                              'showgrid': False,
                              'showline': False,
                              'zeroline': False,
                              'showticklabels': False,
                              'ticks': ''})

    # change layout design and add title
    fig.update_layout(template='seaborn', paper_bgcolor='rgba(0,0,0,0)',
                      plot_bgcolor='rgba(0,0,0,0)',
                      title='Dendogram and allelic differences heatmap')

    # create HTML file with dendogram and heatmap
    plot(fig, filename=output_file, auto_open=False)

    return True


def write_list(input_list, output_file, mode='w'):
    r"""Write list to file.

    Each nested list is joined with '\t' and written as a single line.

    Parameters
    ----------
    input_list : list
        List to write.
    output_file : str
        Path to the output file.
    mode : str
        Write mode ('w' for write, 'a' for append).
    """
    joined_lines = ['\t'.join(l) for l in input_list]
    joined_text = '\n'.join(joined_lines)
    with open(output_file, mode) as outfile:
        outfile.write(joined_text+'\n')


def create_directory(directory_path):
    """Create a directory if it does not exist."""
    if os.path.isdir(directory_path) is False:
        os.mkdir(directory_path)


def divide_results_per_cluster(clusters, allelecall_df,
                               output_directory, concat_file):
    """Divide allele call results into one file per cluster.

    Parameters
    ----------
    clusters : dict
        Dictionary with cluster identifiers as keys and
        a list with the sample identifiers in each cluster
        as values.
    allelecall_df : pandas.core.frame.DataFrame
        Pandas dataframe that represents an allele call matrix.
    output_directory : str
        Path to the output directory to which the files will
        be written.
    concat_file : str
        Path to the file that will contain the concatenated
        and ordered results.

    Returns
    -------
    clusters_files : dict
        Dictionary with cluster identifiers as keys and
        a list that contains the path to the cluster directory
        and the TSV file with allele calling results as values.
    """
    clusters_files = {}
    include_header = True
    for k, v in clusters.items():
        cluster_results = allelecall_df.loc[v]
        # Create subdir to store cluster results
        cluster_dir = os.path.join(output_directory, k)
        create_directory(cluster_dir)

        # save cluster df to cluster dir
        cluster_tsv = os.path.join(cluster_dir, '{0}_allelecall_results.tsv'.format(k))
        cluster_results.to_csv(cluster_tsv, sep='\t')

        clusters_files[k] = [cluster_dir, cluster_tsv]

        # appends cluster profiles to file with all results
        cluster_results.to_csv(concat_file, sep='\t',
                               mode='a', header=include_header)

        include_header = False

    return clusters_files


def create_dm(input_file, output_directory):
    """Compute distance matrix from allele calling results.

    Parameters
    ----------
    input_file : str
        Path to a TSV file with allele call results.
    output_directory : str
        Path to the output directory.

    Returns
    -------
    List with path to a TSV file with the distance matrix
    based on the allelic differences and path to a TSV
    file with the matrix for the number fo shared loci.
    """
    input_basename = os.path.basename(input_file)
    # remove extension that is after last '.'
    input_basename = '.'.join(input_basename.split('.')[0:-1])

    tmp_directory = os.path.join(output_directory, 'tmp')
    create_directory(tmp_directory)

    genome_ids = dm.get_sample_ids(input_file, delimiter='\t')

    np_matrix = dm.tsv_to_nparray(input_file)
    # compute distances per row to avoid heavy memory usage
    rows_indexes = [i for i in range(len(np_matrix))]
    results = dm.compute_distances(rows_indexes, np_matrix,
                                   genome_ids, tmp_directory)

    print('\nCreating distance matrix...', end='')
    # create files with headers
    col_ids = ['FILE'] + genome_ids
    output_pairwise = os.path.join(output_directory,
                                   '{0}_allelic_differences.tsv'.format(input_basename))
    output_p = os.path.join(output_directory,
                            '{0}_shared_loci.tsv'.format(input_basename))

    # import arrays per genome and save to matrix file
    results2 = dm.write_matrices(results, genome_ids,
                                 output_pairwise, output_p, col_ids)

    # add 1 to include header
    symmetric_allelic_differences = dm.symmetrify_matrix(output_pairwise,
                                                         len(genome_ids)+1,
                                                         tmp_directory)
    symmetric_shared_loci = dm.symmetrify_matrix(output_p,
                                                 len(genome_ids)+1,
                                                 tmp_directory)

    shutil.rmtree(tmp_directory)
    
    return [symmetric_allelic_differences, symmetric_shared_loci]


def inter_cluster_stats(stats, distance_matrix):
    """Compute intercluster distance statistics for all clusters.

    Parameters
    ----------
    stats : dict
        Dictionary with the distance statsitics per cluster.
    distance_matrix : str
        Path to a TSV file that contains a distance matrix
        with the number of allelic differences.

    Returns
    -------
    stats : dict
        Dictionary with statistics updated with the
        intercluster distance values.
    """
    cluster_ids = [k for k in stats if 'centroid' not in k]
    current_cols = ['FILE'] + cluster_ids

    stats_df = pd.read_csv(distance_matrix,
                           usecols=current_cols,
                           sep='\t',
                           index_col='FILE')

    outdf = stats_df[~stats_df.index.isin(cluster_ids)]
    # do not compute values if there is only one cluster
    if len(outdf) != 0:
        for e in cluster_ids:
            strain_stats = intercluster_stats(outdf, e)
            stats[e]['min_out'] = strain_stats[0]
            stats[e]['max_out'] = strain_stats[1]
            stats[e]['mean_out'] = strain_stats[2]

    return stats


def boxplot_html(data, output_file, title, xaxis_title,
                 yaxis_title, legend_title, boxplot_min,
                 fontsize):
    """Create HTML file with boxplots.

    Parameters
    ----------
    data : list
        List with boxplot traces.
    output_file : str
        Path to the output HTML file.
    title : str
        Title to display at the top of the page.
    xaxis_title : str
        Title for the xaxis.
    yaxis_title : str
        Title for the yaxis.
    legend_title : str
        Title for the legend.
    """
    # Create HTML with Boxplots
    fig = go.Figure()
    for t in data:
        # only add traces for clusters with more than minimum_n samples
        if boxplot_min is None:
            fig.add_trace(t)
        else:
            if len(t.y) >= int(((boxplot_min*boxplot_min)-boxplot_min)/2):
                fig.add_trace(t)

    fig.update_layout(template='seaborn',
                      title=title,
                      xaxis_title=xaxis_title,
                      yaxis_title=yaxis_title,
                      legend_title=legend_title,
                      font=dict(size=fontsize))

    plot(fig, filename=output_file, auto_open=False)


def heatmap_html(data, output_file, title):
    """Create HTML file with a heatmap.

    Parameters
    ----------
    data : go.Heatmap
        Heatmap trace.
    output_file : str
        Path to the output HTML file.
    title : str
        Title to display at the top of the page.
    """
    # Heatmap
    fig = go.Figure()
    fig.add_trace(data)
    fig.update_layout(template='seaborn', title=title)

    plot(fig, filename=output_file, auto_open=False)


def write_intracluster_stats(data, output_file, headers):
    """Write TSV file with intracluster statistics.

    Parameters
    ----------
    data : dict
        Dictionary with sample identifiers as keys and
        a dictionary with intracluster statistics as values.
    output_file : str
        Path to the output file.
    headers : list
        File headers.
    """
    lines = [headers]
    for h, j in data.items():
        if 'centroid' not in h:
            current_id = h
            min_distance = j['min']
            min_id = min_distance[0][0]
            min_distance = str(min_distance[1])

            max_distance = j['max']
            max_id = max_distance[0][0]
            max_distance = str(max_distance[1])

            mean_distance = j['mean']
            if mean_distance != '':
                mean_distance = str(round(mean_distance, 2))

            lines.append([h, min_distance, min_id,
                          max_distance, max_id,
                          mean_distance])

    write_list(lines, output_file)


def global_stats_lines(data, cluster_id):
    """Write TSV file with cluster statistics for the analysis of all clusters.

    Parameters
    ----------
    data : dict
        Dictionary with sample identifiers as keys and
        a dictionary with intracluster and intercluster
        statistics as values.
    """
    lines = []
    for h, j in data.items():
        if 'centroid' not in h:
            current_id = h
            min_distance = j['min']
            min_id = min_distance[0][0]
            min_distance = str(min_distance[1])

            max_distance = j['max']
            max_id = max_distance[0][0]
            max_distance = str(max_distance[1])

            mean_distance = j['mean']
            if mean_distance != '':
                mean_distance = str(round(mean_distance, 2))

            min_distance_out = j['min_out']
            min_id_out = min_distance_out[0][0]
            min_distance_out = str(min_distance_out[1])

            max_distance_out = j['max_out']
            max_id_out = max_distance_out[0][0]
            max_distance_out = str(max_distance_out[1])

            mean_distance_out = str(round(j['mean_out'], 2))

            lines.append([current_id, min_distance, min_id,
                          max_distance, max_id,
                          mean_distance, min_distance_out,
                          min_id_out, max_distance_out,
                          max_id_out, mean_distance_out,
                          cluster_id])

    return lines


def custom_cgMLST(allelecall_results, core_genome, output_file):
    """Read allele call results for a set of loci.

    Parameters
    ----------
    allelecall_results : str
        Path to a TSV file with allele call results.
    core_genome : list
        Path to file with the list of loci that
        constitute the core genome.
    output_file : str
        Path to the output file.
    """
    # read list of loci in cgMLST
    with open(core_genome, 'r') as infile:
        cgMLST_loci = infile.read().splitlines()

    # read allelecall results for those loci
    allelecall_df = pd.read_csv(allelecall_results,
                                usecols=['FILE']+cgMLST_loci,
                                sep='\t',
                                index_col='FILE')

    # save cgMLST results
    allelecall_df.to_csv(output_file, sep='\t')


def main(input_data, output_directory, clusters, cpu_cores,
         dendogram_threshold, dataset_name, core_genome, minimum_n,
         plot_fontsize):

    # Create output directory
    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    clusters_dict = {}
    if clusters is not None:
        # read clusters
        clusters_lines = read_tabular(clusters)
        for l in clusters_lines:
            clusters_dict.setdefault(l[1], []).append(l[0])

    # sort clusters by decreasing number of elements
    clusters_size = [(k, len(v))
                     for k, v in clusters_dict.items()]
    clusters_size = sorted(clusters_size,
                           key=lambda x: x[1], reverse=True)
    ordered_clusters = {c[0]: clusters_dict[c[0]]
                        for c in clusters_size}

    # get input basename to use as file prefix
    input_basename = os.path.basename(input_data)
    input_basename = input_basename.split('.tsv')[0]

    # mask input matrix
    print('Masking input matrix...', end='')
    masked_results = os.path.join(output_directory, input_basename+'_masked.tsv')
    nostdout(mm.main, [input_data, masked_results, None])
    print('done')

    # import masked matrix
    allelecall_df = pd.read_csv(masked_results, sep='\t', index_col='FILE')

    all_results_dir = os.path.join(output_directory, 'all_results')
    create_directory(all_results_dir)
    all_results_file = os.path.join(all_results_dir, 'all_results.tsv')
    if clusters is not None:
        clusters_dirs = divide_results_per_cluster(ordered_clusters, allelecall_df,
                                                   output_directory, all_results_file)
        clusters_dirs['all_results'] = [all_results_dir, all_results_file]

    # determine cgMLST for all clusters
    for k, v in clusters_dirs.items():
        print('Determining cgMLST for {0}...'.format(k), end='')
        cg_dir = os.path.join(v[0], '{0}_cgMLST'.format(k))
        cg_file = os.path.join(cg_dir, 'cgMLST.tsv')
        # determine cgMLST for each cluster
        if core_genome is None:
            current_cg = nostdout(cg.main, [v[1], cg_dir, 1.0, [], []])
        # use cgMLST provided by user
        else:
            create_directory(cg_dir)
            custom_cgMLST(v[1], core_genome, cg_file)
        print('done.')

        clusters_dirs[k].append(cg_file)

    # create directories to store files with distances
    for k, data in clusters_dirs.items():
        # wgMLST
        wgMLST_distance_dir = os.path.join(data[0], 'wgMLST_{0}_distances'.format(k))
        create_directory(wgMLST_distance_dir)

        # cgMLST
        cgMLST_distance_dir = os.path.join(data[0], 'cgMLST_{0}_distances'.format(k))
        create_directory(cgMLST_distance_dir)

        clusters_dirs[k].extend([wgMLST_distance_dir, cgMLST_distance_dir])

    # running the main function from the distance_matrix.py script
    # originates a memory leak. Could not determine the cause.
    # running other functions that are not the main function does
    # not cause a memory leak

    # determine distance matrices
    processed = 0
    for k, v in clusters_dirs.items():
        print('Computing distance matrix for {0}...'.format(k), end='')
        # wgMLST
        wgMLST_res = nostdout(create_dm, [v[1], v[3]])
        # cgMLST
        cgMLST_res = nostdout(create_dm, [v[2], v[4]])

        clusters_dirs[k].extend([wgMLST_res[0], cgMLST_res[0]])
        print('done.')

    # determine intra-cluster stats for all clusters
    all_stats = {}
    for k, v in clusters_dirs.items():
        print('Computing distance statistics for {0}...'.format(k), end='')
        # determine wg and cg stats for all samples
        if k == 'all_results':
            current_clusters = ordered_clusters
        # determine wg and cg stats per cluster
        else:
            current_clusters = {k: ordered_clusters[k]}

        # wgMLST
        wgMLST_file = v[5]
        wg_stats = clusters_stats(wgMLST_file, current_clusters)

        # cgMLST
        cgMLST_file = v[6]
        cg_stats = clusters_stats(cgMLST_file, current_clusters)

        all_stats[k] = [wg_stats, cg_stats]
        print('done.')

    # create plots for all clusters
    traces = {}
    for k, v in clusters_dirs.items():

        if k == 'all_results':
            current_clusters = ordered_clusters
        else:
            current_clusters = {k: ordered_clusters[k]}

        # wgMLST
        wgMLST_file = v[5]
        wg_boxplot_trace = clusters_boxplot(wgMLST_file, current_clusters)
        wg_heatmap_trace = clusters_heatmap(wgMLST_file)

        # cgMLST
        cgMLST_file = v[6]
        cg_boxplot_trace = clusters_boxplot(cgMLST_file, current_clusters)
        cg_heatmap_trace = clusters_heatmap(cgMLST_file)

        traces[k] = [[wg_boxplot_trace, wg_heatmap_trace],
                     [cg_boxplot_trace, cg_heatmap_trace]]

    # determine distances between cluster centroids
    # get centroids
    # wgMLST
    wgMLST_centroids = {v['centroids']: k
                        for k, v in all_stats['all_results'][0].items()}

    # only compute if there is more than one cluster
    if len(wgMLST_centroids) > 1:
        distance_matrix = clusters_dirs['all_results'][5]
        for k, v in wgMLST_centroids.items():
            centroids_ids = list(wgMLST_centroids.keys())
            centroids_ids.remove(k)
            centroid_stats = centroids_inter_dists(distance_matrix, centroids_ids, k)
            all_stats['all_results'][0][v]['centroid_min'] = centroid_stats[0]
            all_stats['all_results'][0][v]['centroid_max'] = centroid_stats[1]
            all_stats['all_results'][0][v]['centroid_mean'] = centroid_stats[2]

    # cgMLST
    cgMLST_centroids = {v['centroids']: k
                        for k, v in all_stats['all_results'][1].items()}

    # only compute if there is more than one cluster
    if len(cgMLST_centroids) > 1:
        distance_matrix = clusters_dirs['all_results'][6]
        for k, v in cgMLST_centroids.items():
            centroids_ids = list(cgMLST_centroids.keys())
            centroids_ids.remove(k)
            centroid_stats = centroids_inter_dists(distance_matrix, centroids_ids, k)
            all_stats['all_results'][1][v]['centroid_min'] = centroid_stats[0]
            all_stats['all_results'][1][v]['centroid_max'] = centroid_stats[1]
            all_stats['all_results'][1][v]['centroid_mean'] = centroid_stats[2]

    # get inter-cluster values for all samples
    wgMLST_file = clusters_dirs['all_results'][5]
    cgMLST_file = clusters_dirs['all_results'][6]
    clusters = [k for k in all_stats if k != 'all_results']
    for c in clusters:
        # wgMLST
        data = all_stats['all_results'][0][c]
        res = inter_cluster_stats(data, wgMLST_file)
        all_stats['all_results'][0][c] = res
        # cgMLST
        data = all_stats['all_results'][1][c]
        res = inter_cluster_stats(data, cgMLST_file)
        all_stats['all_results'][1][c] = res

    # Create HTML files with plots
    print('Creating TSV files with distance statistics and HTML files...', end='')
    boxplot_title = '{0} comparison at {1}MLST level (number of allelic differences)'
    heatmap_title = 'Distance matrix at {0}MLST level (number of allelic differences)'
    for k, v in traces.items():
        if k == 'all_results':
            boxplot_min = minimum_n
        else:
            boxplot_min = None

        outdir = os.path.join(output_directory, k)
        # boxplots
        # wgMLST
        output_html = os.path.join(outdir,
                                   'wgMLST_clusters_boxplots.html')
        boxplot_html(v[0][0], output_html, boxplot_title.format(dataset_name, 'wg'),
                     dataset_name, 'Number of allelic differences',
                     dataset_name, boxplot_min, plot_fontsize)
        # cgMLST
        output_html = os.path.join(outdir,
                                   'cgMLST_clusters_boxplots.html')
        boxplot_html(v[1][0], output_html, boxplot_title.format(dataset_name, 'cg'),
                     dataset_name, 'Number of allelic differences',
                     dataset_name, boxplot_min, plot_fontsize)

        # heatmaps
        # wgMLST
        output_html = os.path.join(outdir,
                                   'wgMLST_heatmap.html')
        heatmap_html(v[0][1], output_html, heatmap_title.format('wg'))
        # cgMLST
        output_html = os.path.join(outdir,
                                   'cgMLST_heatmap.html')
        heatmap_html(v[1][1], output_html, heatmap_title.format('cg'))

    # Create plot with boxplots for the separate analysis of each cluster
    # this is different than the boxplots for the analysis of the complete results
    # these boxplots cannot be directly compared because the set of loci in each might be different
    wgMLST_fig = go.Figure()
    wgMLST_html = os.path.join(output_directory, 'wgMLST_clusters_boxplots.html')
    cgMLST_fig = go.Figure()
    cgMLST_html = os.path.join(output_directory, 'cgMLST_clusters_boxplots.html')

    wgMLST_traces = []
    cgMLST_traces = []
    for k, v in traces.items():
        if k != 'all_results':
            wgtrace = v[0][0][0]
            cgtrace = v[1][0][0]
            if len(wgtrace.y) >= int(((minimum_n*minimum_n)-minimum_n)/2):
                wgMLST_traces.append(wgtrace)
                cgMLST_traces.append(cgtrace)

    for i in range(len(wgMLST_traces)):
        wgMLST_traces[i]['marker']['color'] = '#2b8cbe'
        cgMLST_traces[i]['marker']['color'] = '#2b8cbe'
        if wgMLST_traces[i]['fillcolor'] != 'rgba(0,0,0,0)':
            wgMLST_traces[i]['fillcolor'] = '#2b8cbe'
            cgMLST_traces[i]['fillcolor'] = '#2b8cbe'

    for t in wgMLST_traces:
        wgMLST_fig.add_trace(t)
    wgMLST_fig_title = '{0} comparison at wgMLST level (number of allelic differences)'.format(dataset_name)
    wgMLST_fig.update_layout(template='seaborn',
                             title=wgMLST_fig_title,
                             xaxis_title=dataset_name,
                             yaxis_title='Number of allelic differences',
                             legend_title=dataset_name,
                             font=dict(size=plot_fontsize))
    plot(wgMLST_fig, filename=wgMLST_html, auto_open=False)

    for t in cgMLST_traces:
        cgMLST_fig.add_trace(t)
    cgMLST_fig_title = '{0} comparison at cgMLST level (number of allelic differences)'.format(dataset_name)
    cgMLST_fig.update_layout(template='seaborn',
                             title=cgMLST_fig_title,
                             xaxis_title=dataset_name,
                             yaxis_title='Number of allelic differences',
                             legend_title=dataset_name,
                             font=dict(size=plot_fontsize))
    plot(cgMLST_fig, filename=cgMLST_html, auto_open=False)

    # Create output files with stats
    file_headers = ['FILE', 'min_distance', 'min_ID', 'max_distance',
                    'max_ID', 'mean_distance', 'min_out_distance',
                    'min_out_ID', 'max_out_distance', 'max_out_ID',
                    'mean_out_distance', 'cluster']
    cluster_headers = file_headers[0:6]
    for k, v in all_stats.items():
        if k != 'all_results':
            # wgMLST
            outfile = os.path.join(output_directory, k, 'wgMLST_stats.tsv')
            write_intracluster_stats(v[0][k], outfile, cluster_headers)
            
            # cgMLST
            outfile = os.path.join(output_directory, k, 'cgMLST_stats.tsv')
            write_intracluster_stats(v[1][k], outfile, cluster_headers)

    # Create output files with global stats
    global_stats = all_stats['all_results']
    clusters_ids = list(global_stats[0].keys())
    wgMLST_lines = [file_headers]
    cgMLST_lines = [file_headers]
    for k in clusters_ids:
        # wgMLST
        wgMLST_stats = global_stats[0][k]
        current_lines = global_stats_lines(wgMLST_stats, k)
        wgMLST_lines.extend(current_lines)

        # cgMLST
        cgMLST_stats = global_stats[1][k]
        current_lines = global_stats_lines(cgMLST_stats, k)
        cgMLST_lines.extend(current_lines)

    outfile = os.path.join(output_directory, 'all_results', 'wgMLST_stats.tsv')
    write_list(wgMLST_lines, outfile)
    outfile = os.path.join(output_directory, 'all_results', 'cgMLST_stats.tsv')
    write_list(cgMLST_lines, outfile)

    # define linkage function to use
    # Perform hierarchical/agglomerative clustering
    # single/min/nearest linkage on the condensed distance matrix
    linkage_function = lambda x: linkage(x, 'single')

    # distfun simply reads the distance matrix that was computed
    # the implementation of ff.create_dendrogram applied the distfun
    # to the whole array. By defining a distfun that returns the passed
    # argument we can pass a precomputed distance matrix
    distance_function = simply_return

    for k, v in clusters_dirs.items():
        # wgMLST
        distance_matrix_file = v[5]
        output_file = os.path.join(v[0], 'wgMLST_dendogram.html')
        res = cluster_dendogram(distance_matrix_file, output_file,
                                linkage_function, distance_function)

        # cgMLST
        distance_matrix_file = v[6]
        output_file = os.path.join(v[0], 'cgMLST_dendogram.html')
        res = cluster_dendogram(distance_matrix_file, output_file,
                                linkage_function, distance_function)

    print('done.')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-data', type=str,
                        required=True, dest='input_data',
                        help='Path to a TSV file that contains '
                             'allele call results.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('-c', '--clusters', type=str,
                        required=True, dest='clusters',
                        help='Path to a TSV file with sample '
                             'identifiers and the cluster the samples '
                             'belong to (sample ids in the first column '
                             'and cluster ids in the second column). '
                             'The analysis also groups all samples '
                             'into a global cluster that will be analysed '
                             'in the same way as predefined clusters.')

    parser.add_argument('-cpu', '--cpu-cores', type=int,
                        required=False, default=1,
                        dest='cpu_cores',
                        help='Number of CPU cores to use.')

    # currently unused
    parser.add_argument('-dt', '--dendogram-threshold', type=float,
                        required=False, default=None,
                        dest='dendogram_threshold',
                        help='Threshold to define dendogram groups.')

    parser.add_argument('-dn', '--dataset-name', type=str,
                        required=False, default='My dataset',
                        dest='dataset_name',
                        help='Name of the dataset.')

    parser.add_argument('-cg', '--core-genome', type=str,
                        required=False,
                        dest='core_genome',
                        help='Provide list of loci in the cgMLST.'
                             'Path to TXT file with list of loci '
                             'that constitute the core genome '
                             '(one locus per line). The analysis '
                             'will use the provided cgMLST for all '
                             'clusters. It determines the cgMLST '
                             'for each cluster when no cgMLST is '
                             'provided.')

    parser.add_argument('-mn', '--minimum-n', type=int,
                        required=False, default=0,
                        dest='minimum_n',
                        help='Minimum number of isolates that a cluster '
                             'must include to create a boxplot for that '
                             'cluster.')

    parser.add_argument('-pf', '--plot-fontsize', type=int,
                        required=False, default=14,
                        dest='plot_fontsize',
                        help='Fontsize for text in HTML files.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
