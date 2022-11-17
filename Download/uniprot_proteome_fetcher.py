#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a TSV file downloaded from UniProt with
information about a set of proteomes and downloads the
proteomes listed in the file.

Code documentation
------------------
"""


import os
import sys
import csv
import time
import socket
import argparse
import urllib.request
import concurrent.futures
from itertools import repeat


# set socket timeout for urllib calls
socket.setdefaulttimeout(30)

# URL template for proteome download
proteome_template_url = 'https://www.uniprot.org/uniprot/?query=proteome:{0}&format=fasta&compress=yes'


def download_file(url, file_name, retry):
    """Accept a URL to download a file.

    Parameters
    ----------
    url : str
        An url to download a file.
    file_name : str
        The name of the file to be downloaded.
    retry : int
        Maximum number of retries if download fails.

    Returns
    -------
    response : str
        A string indicating that the download failed or
        an object with the response information for a
        successful download.
    """

    tries = 0
    while tries < retry:
        try:
            response = urllib.request.urlretrieve(url, file_name)
            break
        except Exception:
            response = 'Failed: {0}'.format(file_name)
            tries += 1
            print('Retrying {0} ...{1}'.format(file_name.split('/')[-1], tries))
            time.sleep(1)

    return response


def main(input_table, output_directory, threads, retry):

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # open table downloaded from Uniprot
    with open(input_table, 'r') as table:
        reader = csv.reader(table, delimiter=',')
        # exclude header
        next(reader, None)
        proteomes = [row[0].split("\t")[0] for row in reader]

    # Build the Uniprot URLs
    urls = [proteome_template_url.format(proteome_id)
            for proteome_id in proteomes]

    files_number = len(urls)
    if files_number == 0:
        sys.exit('No valid URLs.')

    local_filenames = ['{0}.fasta.gz'.format(proteome_id)
                       for proteome_id in proteomes]
    local_filepaths = [os.path.join(output_directory, filename)
                       for filename in local_filenames]

    print('\nStarting download of {0} proteomes...'.format(len(urls)))
    start = time.time()

    # We can use a with statement to ensure threads are cleaned up promptly
    failures = []
    success = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_file, urls, local_filepaths, repeat(retry)):
            if 'Failed' in res:
                failures.append(res)
            else:
                success += 1
                print('\r', 'Downloaded {0}/{1}'.format(success, files_number), end='')

    print('\nFailed download for {0} files.'.format(len(failures)))

    end = time.time()
    delta = end - start
    minutes = int(delta/60)
    seconds = delta % 60
    print('\nFinished downloading {0}/{1} fasta.gz files.'
          '\nElapsed Time: {2}m{3:.0f}s'
          ''.format(files_number-len(failures),
                    files_number, minutes, seconds))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--input_table', type=str,
                        required=True, dest='input_table',
                        help='TSV file downloaded from UniProt '
                             'that contains list of proteomes.')

    parser.add_argument('-o', '--output_directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory where downloaded '
                             'files will be stored.')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for concurrent download.')

    parser.add_argument('-r', '--retry', type=int,
                        required=False, dest='retry',
                        default=7,
                        help='Maximum number of retries when a '
                             'download fails.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
