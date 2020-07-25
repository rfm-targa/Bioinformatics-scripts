#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

    This script accepts a Genome Assembly and Annotation report
    table from the NCBI and downloads the genomes/assemblies of
    the samples listed in the table. If the genomes/assemblies
    are not in the refseq FTP, it tries to download from Genbank.
    After failing a download, the script keeps trying up to 7 times.
"""


import os
import csv
import time
import socket
import argparse
import urllib.request
import concurrent.futures


socket.setdefaulttimeout(30)


def download_assembly(url, file_name):
    """ Accepts a url to download a file, retrying up to 7 times
        if previous attempts are not successful.

        Args:
            url (str): an url to download a file.
            file_name (str): the identifier of the file to be downloaded.

        Returns:
            response: a string indicating that the download failed or
            an object with the response information for the successful
            download.
    """

    tries = 0
    download_tries = 7
    while tries < download_tries:
        try:
            response = urllib.request.urlretrieve(url, file_name)
            tries = 7
        except Exception:
            response = 'Failed: {0}'.format(file_name)
            tries += 1
            print('Retrying {0} ...{1}'.format(file_name.split('/')[-1], tries))

    return response


def main(input_table, output_directory, file_extension):

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # open table downloaded from NCBI
    with open(input_table, 'r') as table:
        lines = list(csv.reader(table, delimiter=','))

    # get urls for samples that have refseq ftp path
    refseq_urls = [line[15] for line in lines[1:] if line[15].strip() != '']
    refseq_assemblies_ids = [url.split('/')[-1] for url in refseq_urls]

    # try to get genbank urls for samples that had no refseq ftp path
    genbank_urls = [line[14] for line in lines[1:] if line[15].strip() == '' and line[14].strip() != '']
    genbank_assemblies_ids = [url.split('/')[-1] for url in genbank_urls]

    urls = refseq_urls + genbank_urls
    assemblies_ids = refseq_assemblies_ids + genbank_assemblies_ids

    ftp_urls = []
    for url in range(len(urls)):
        ftp_url = '{0}/{1}_{2}'.format(urls[url],
            assemblies_ids[url], file_extension)
        ftp_urls.append(ftp_url)

    files_number = len(ftp_urls)

    assemblies_ids = ['{0}.fasta.gz'.format(url) for url in assemblies_ids]
    assemblies_ids = [os.path.join(output_directory, file_name)
                      for file_name in assemblies_ids]

    print('\nStarting download of {0} fasta.gz files...'.format(len(ftp_urls)))
    start = time.time()

    # We can use a with statement to ensure threads are cleaned up promptly
    failures = []
    success = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=15) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_assembly, ftp_urls, assemblies_ids):
            failures.append(res)
            success += 1
            print('\r', 'Downloaded {0}'.format(success), end='')

    failures_number = len([res for res in failures if 'Failed' in res])

    end = time.time()
    delta = end - start
    minutes = int(delta/60)
    seconds = delta % 60
    print('\nFinished downloading {0}/{1} fasta.gz files.'
          '\nElapsed Time: {2}m{3:.0f}s'
          ''.format(files_number-failures_number,
                    files_number, minutes, seconds))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-t', '--input_table', type=str,
                        required=True, dest='input_table',
                        help='TSV file downloaded from the NCBI '
                             'Genome Assembly and Annotation report.')

    parser.add_argument('-o', '--output_directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory where downloaded '
                             'files will be stored.')

    parser.add_argument('--fe', '--file_extension', type=str,
                        required=False, default='genomic.fna.gz',
                        dest='file_extension',
                        help='')

    args = parser.parse_args()

    return [args.input_table, args.output_directory, args.file_extension]


if __name__ == '__main__':

    args = parse_arguments()
    main(args[0], args[1], args[2])
