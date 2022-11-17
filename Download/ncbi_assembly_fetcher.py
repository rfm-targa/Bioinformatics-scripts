#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a Genome Assembly and Annotation report
table from the NCBI and downloads the genome assemblies of
the samples listed in the table.

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


# url to download assembly_summary_refseq.txt
assembly_summary_refseq = 'https://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt'


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


def main(input_table, output_directory, file_extension,
         ftp, threads, organism_name, retry):

    organism_name = organism_name.lower()

    if not os.path.isdir(output_directory):
        os.mkdir(output_directory)

    # open table downloaded from NCBI
    with open(input_table, 'r') as table:
        lines = list(csv.reader(table, delimiter=','))

    # get urls for samples that have refseq ftp path
    refseq_urls = []
    refseq_assemblies_ids = []
    if 'refseq' in ftp:
        refseq_urls = [line[15]
                       for line in lines[1:]
                       if line[15].strip() != '']

        refseq_assemblies_ids = [url.split('/')[-1]
                                 for url in refseq_urls]

        print('Found URLs for {0} strains in RefSeq.'
              ''.format(len(refseq_urls)))

    genbank_urls = []
    genbank_assemblies_ids = []
    if 'genbank' in ftp:
        # get genbank urls
        # only get if sample is not in refseq list
        genbank_urls = [line[14]
                        for line in lines[1:]
                        if line[14].strip() != '' and line[15] not in refseq_urls]

        genbank_assemblies_ids = [url.split('/')[-1]
                                  for url in genbank_urls]

        print('Found URLs for {0} strains in GenBank.'
              ''.format(len(genbank_urls)))

    print('Downloading assembly_summary_refseq...', end='')
    assembly_summary_refseq_local = os.path.join(output_directory,
                                                 'assembly_summary_refseq.txt')
    download_file(assembly_summary_refseq, assembly_summary_refseq_local, retry)
    print('done.')

    # open assembly_summary_refseq table 
    with open(assembly_summary_refseq_local, 'r') as table1:
        assembly_summary_lines = list(csv.reader(table1, delimiter='\t'))

    os.remove(assembly_summary_refseq_local)

    # processing assembly_summary_refseq table
    # filtering for the organism_name given as argument
    assembly_latest = []
    for row in assembly_summary_lines[1:]:
        if organism_name in row[7].lower():
            assembly_latest.append(row[0])

    # removing suppressed genomes from refseq_assemblies_ids and refseq_urls
    suppressed = []
    for i in refseq_assemblies_ids:
        processed_id = '_'.join(i.split('_')[0:2])

        if processed_id not in assembly_latest:
            suppressed.append(refseq_assemblies_ids.index(i))

    print('Found {0} ids suppressed from RefSeq.'.format(len(suppressed)))

    # remove suppressed from lists
    refseq_assemblies_ids = [j
                             for i, j in enumerate(refseq_assemblies_ids)
                             if i not in suppressed]
    refseq_urls = [j
                   for i, j in enumerate(refseq_urls)
                   if i not in suppressed]

    urls = refseq_urls + genbank_urls
    assemblies_ids = refseq_assemblies_ids + genbank_assemblies_ids

    # construct FTP URLs
    ftp_urls = []
    for i, url in enumerate(urls):
        ftp_url = '{0}/{1}_{2}'.format(url,
                                       assemblies_ids[i],
                                       file_extension)
        ftp_urls.append(ftp_url)

    files_number = len(ftp_urls)
    if files_number == 0:
        sys.exit('No valid ftp links after scanning table.')

    assemblies_ids = ['{0}/{1}_{2}'.format(output_directory, url, file_extension)
                      for url in assemblies_ids]

    print('\nStarting download of {0} {1} files...'
          ''.format(files_number, file_extension))

    start = time.time()
    # We can use a with statement to ensure threads are cleaned up promptly
    failures = []
    success = 0
    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        # Start the load operations and mark each future with its URL
        for res in executor.map(download_file, ftp_urls, assemblies_ids, repeat(retry)):
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

    parser.add_argument('-t', '--input-table', type=str,
                        required=True, dest='input_table',
                        help='TSV file downloaded from the '
                             'NCBI Genome Assembly and Annotation '
                             'report.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the directory to which '
                             'downloaded files will be stored.')

    parser.add_argument('-n', '--organism-name', type=str,
                        required=True, dest='organism_name',
                        help='Organism name. For example: '
                             'Streptococcus pneumoniae.')

    parser.add_argument('--fe', '--file-extension', type=str,
                        required=False,
                        choices=['genomic.fna.gz', 'assembly_report.txt',
                                 'assembly_status.txt', 'cds_from_genomic.fna.gz',
                                 'feature_count.txt.gz', 'feature_table.txt.gz',
                                 'genomic.gbff.gz', 'genomic.gff.gz',
                                 'genomic.gtf.gz', 'protein.faa.gz',
                                 'protein.gpff.gz', 'rna_from_genomic.fna.gz',
                                 'translated_cds.faa.gz'],
                        default='genomic.fna.gz',
                        dest='file_extension',
                        help='Choose file type to download through '
                             'extension.')

    parser.add_argument('--ftp', type=str, required=False,
                        choices=['refseq+genbank', 'refseq', 'genbank'],
                        default='refseq+genbank', dest='ftp',
                        help='The script can search for the files to '
                             'download in RefSeq or Genbank or both '
                             '(will only search in Genbank if download '
                             'from RefSeq fails).')

    parser.add_argument('-th', '--threads', type=int,
                        required=False, default=2,
                        dest='threads',
                        help='Number of threads for download.')

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
