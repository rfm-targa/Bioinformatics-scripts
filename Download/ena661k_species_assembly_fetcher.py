#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script downloads genome assemblies from the
ENA661k project (https://doi.org/10.1101/2021.03.02.433662).
It can be used to download all genome assemblies for
a species and the assemblies to download can be filtered
based on a set of quality criteria.

"""


import os
import sys
import csv
import time
import argparse
import urllib.request
from multiprocessing import TimeoutError
from multiprocessing.pool import ThreadPool


# change csv module settings to read huge files
maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)


# ftp path
ebi_ftp = 'http://ftp.ebi.ac.uk'


def read_table(file_path, delimiter='\t'):
    """ Reads a tabular file.

        Parameters
        ----------
        file_path : str
            Path to the tabular file.
        delimiter : str
            Field delimiter.

        Returns
        -------
        lines : list
            List that contains one sublist per
            line in the input tabular file.
    """

    with open(file_path, 'r') as infile:
        lines = list(csv.reader(infile, delimiter=delimiter))

    return lines


def input_timeout(func, input_args, timeout):
    """ Adds timeout feature to a function call.

        Parameters
        ----------
        func
            Function to call.
        input_args : list
            List with the ordered arguments to pass
            to the function.
        timeout : int
            Maximum number of seconds that the process will
            wait for a result.

        Returns
        -------

        Raises
        ------
        SystemExit
            - If there is no return value from the function
              call in the defined timeout.
    """

    pool = ThreadPool(processes=1)
    answer = pool.apply_async(func, args=[*input_args])

    try:
        return answer.get(timeout=timeout)
    except TimeoutError as e:
        sys.exit('Timed out.')


def download_ftp_file(file_url, out_file, minimum_file_size, timeout):
    """ Downloads a file from a FTP server.

        Parameter
        ---------
        file_url : str
            FTP path to the file to download.
        out_file : str
            Local path to the file that will be downloaded.

        Returns
        -------
        downloaded : bool
            True if the file was successfully downloaded,
            False otherwise.
    """

    tries = 0
    downloaded = False
    while tries <= 5 and downloaded is False:
        try:
            res = input_timeout(urllib.request.urlretrieve,
                                [file_url, out_file],
                                timeout=timeout)
            file_size = os.path.getsize(out_file)
            if file_size < minimum_file_size:
                os.remove(out_file)
                raise ValueError
            else:
                downloaded = True
        except:
            if os.path.isfile(out_file) is True:
                os.remove(out_file)
            time.sleep(5)
            print('Retrying {0}'.format(file_url))
        tries += 1

    return downloaded


def download_assemblies(sample_ids, sample_paths, output_directory,
                        minimum_file_size, timeout):
    """ Downloads a set of assemblies from the FTP server of
        the "ENA2018-bacteria-661k" study.

        Parameters
        ----------
        sample_ids : list
            List with the identifiers of the samples/assemblies
            to download.
        sample_paths : dict
            Dictionary with sample/assemblies identifiers as
            keys and FTP paths as values.
        output_directory : str
            Path to the output directory.

        Returns
        -------
        failed : int
            Number of failed downloads.
    """

    # list files in output directory
    local_files = os.listdir(output_directory)

    # create URLs to download
    remote_urls = []
    for sid in sample_ids:
        sample_basename = sample_paths[sid].split('/')[-1]
        # do not download files that have already been downloaded
        if sample_basename not in local_files and sample_basename.split('.gz')[0] not in local_files:
            sample_file = os.path.join(output_directory, sample_basename)
            sample_url = ebi_ftp + sample_paths[sid]
            remote_urls.append([sample_url, sample_file])

    if len(remote_urls) < len(sample_ids):
        print('{0} assemblies had already been downloaded.'
              ''.format(len(sample_ids)-len(remote_urls)))

    print('\nDownloading {0} assemblies...'.format(len(remote_urls)))
    failed = 0
    downloaded = 0
    for url in remote_urls:
        res = download_ftp_file(url[0], url[1], minimum_file_size, timeout)
        if res is True:
            downloaded += 1
            print('Downloaded {0} ({1}/{2})'.format(url[1],
                                                    downloaded,
                                                    len(remote_urls)))
        else:
            failed += 1

    return failed


def main(metadata_table, paths_table, species_name, output_directory,
         ftp_download, abundance, genome_size, size_threshold,
         max_contig_number, minimum_file_size, timeout,
         mlst_species, known_st, any_quality):

    # read file with metadata
    metadata_lines = read_table(metadata_table)

    metadata_header = metadata_lines[0]

    # select lines based on species name
    species_index = metadata_header.index('species')
    species_lines = [line
                     for line in metadata_lines[1:]
                     if line[species_index] == species_name]

    print('\nFound {0} samples for species={1}.'
          ''.format(len(species_lines), species_name))

    if len(species_lines) == 0:
        print('Did not find matches for {0}.'.format(species_name))
        print('Please provide a valid species name.')
        sys.exit(0)

    # filter based on genome size
    if genome_size is not None:
        bot_limit = genome_size - (genome_size*size_threshold)
        top_limit = genome_size + (genome_size*size_threshold)
        size_index = metadata_header.index('total_length')
        species_lines = [line
                         for line in species_lines
                         if int(line[size_index]) >= bot_limit
                         and int(line[size_index]) <= top_limit]

        print('{0} with genome size >= {1} and '
              '<= {2}.'.format(len(species_lines), bot_limit, top_limit))

    # filter based on species abundance
    if abundance is not None:
        abundance_index = metadata_header.index('adjust_abundance')
        species_lines = [line
                         for line in species_lines
                         if float(line[abundance_index]) >= abundance]

        print('{0} with abundance >= {1}.'.format(len(species_lines),
                                                  abundance))

    # filter based on number of contigs
    if max_contig_number is not None:
        contigs_index = metadata_header.index('total_contigs')
        species_lines = [line
                         for line in species_lines
                         if int(line[contigs_index]) <= max_contig_number]

        print('{0} with <= {1} contigs.'.format(len(species_lines),
                                                max_contig_number))

    # filter based on MLST species
    if mlst_species is not None:
        st_species_index = metadata_header.index('mlst-species')
        species_lines = [line
                         for line in species_lines
                         if line[st_species_index] == mlst_species]

        print('{0} with MLST species == {1}.'.format(len(species_lines),
                                                     mlst_species))

    # filter based on known ST
    if known_st is True:
        st_index = metadata_header.index('mlst')
        species_lines = [line
                         for line in species_lines
                         if line[st_index] != '-']

        print('{0} with known ST.'.format(len(species_lines)))

    # filter based on quality level
    if any_quality is False:
        quality_index = metadata_header.index('high_quality')
        species_lines = [line
                         for line in species_lines
                         if line[quality_index] == 'TRUE']

        print('{0} with high quality.'.format(len(species_lines)))

    # get sample identifiers
    sample_ids = [line[0] for line in species_lines]
    print('Selected {0} samples/assemblies that meet filtering '
          'criteria.'.format(len(sample_ids)))

    if len(sample_ids) == 0:
        sys.exit('Did not find samples/assemblies that passed '
                 'filtering criteria.')

    # create output directory
    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)

    selected_file = os.path.join(output_directory, 'selected_samples.tsv')
    with open(selected_file, 'w') as outfile:
        selected_lines = ['\t'.join(line)
                          for line in [metadata_header]+species_lines]
        selected_text = '\n'.join(selected_lines)
        outfile.write(selected_text)

    if ftp_download is True:

        # read table with FTP paths
        ftp_lines = read_table(paths_table)
        sample_paths = {l[0]: l[1].split('/ebi/ftp')[1] for l in ftp_lines}

        # download assemblies
        download_assemblies(sample_ids, sample_paths, output_directory,
                            minimum_file_size, timeout)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-m', '--metadata-table', type=str,
                        required=True, dest='metadata_table',
                        help='Summarised QC and general characterisation for '
                             'each assembly in the "File4_QC_characterisati'
                             'on_661K" file from the study "Exploring '
                             'bacterial diversity via a curated and '
                             'searchable snapshot of archived DNA '
                             'sequences" (https://doi.org/10.1101/2021.03'
                             '.02.433662).')

    parser.add_argument('-p', '--paths-table', type=str,
                        required=True, dest='paths_table',
                        help='File with sample identifier to FTP path '
                             'mapping (available at http://ftp.ebi.ac.uk/'
                             'pub/databases/ENA2018-bacteria-661k/).')

    parser.add_argument('-s', '--species-name', type=str,
                        required=True, dest='species_name',
                        help='Name of the species. Must match one of the '
                             'species names in the "species" column in the '
                             'metadata table.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory.')

    parser.add_argument('--ftp-download', action='store_true',
                        required=False, dest='ftp_download',
                        help='If the assemblies from the selected samples'
                             'should be downloaded.')

    parser.add_argument('-a', '--abundance', type=float,
                        required=False,
                        dest='abundance',
                        help='Minimum species abundance. Samples with species'
                        ' abundance below this value are not selected.')

    parser.add_argument('-gs', '--genome-size', type=int,
                        required=False,
                        dest='genome_size',
                        help='Expected genome size.')

    parser.add_argument('-st', '--size-threshold', type=float,
                        required=False,
                        dest='size_threshold',
                        help='Genome size can vary in size +/- this value.')

    parser.add_argument('-mc', '--max-contig-number', type=int,
                        required=False,
                        dest='max_contig_number',
                        help='Maximum number of contigs. Assemblies with '
                             'a number of contigs greater than this value '
                             'are not selected.')

    parser.add_argument('-mfs', '--minimum-file-size', type=int,
                        required=False, default=500000,
                        dest='minimum_file_size',
                        help='Minimum file size for downloaded files. '
                             'The process will retry downloading files '
                             'that are smaller than this value.')

    parser.add_argument('-t', '--timeout', type=int,
                        required=False, default=30,
                        dest='timeout',
                        help='Maximum number of seconds for the '
                             'download of a file.')

    parser.add_argument('--mlst-species', type=str,
                        required=False,
                        dest='mlst_species',
                        help='The species predicted by the MLST tool.')

    parser.add_argument('--known-st', action='store_true',
                        required=False,
                        dest='known_st',
                        help='If the samples must have a known ST.'
                             'Invalid or unkown STs will be "-".')

    parser.add_argument('--any-quality', action='store_true',
                        required=False,
                        dest='any_quality',
                        help='Download all assemblies, even the ones '
                             'that are not high quality.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
