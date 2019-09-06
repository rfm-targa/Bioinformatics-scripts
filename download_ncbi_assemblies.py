#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION

"""

import os
import csv
import time
import socket
import urllib.request
import concurrent.futures

socket.setdefaulttimeout(30)

workdir = os.path.join(os.getcwd(), 'downloaded_assemblies')

# open table downloaded from NCBI
with open('genomes_proks154.txt', 'r') as table:
    lines = list(csv.reader(table, delimiter='\t'))

urls = [line[17] for line in lines[1:] if line[17] != '-']
assemblies_ids = [url.split('/')[-1] for url in urls]

ftp_urls = []
for url in range(len(urls)):
    ftp_url = '{0}/{1}_genomic.fna.gz'.format(urls[url], assemblies_ids[url])
    ftp_urls.append(ftp_url)

assemblies_ids = ['{0}.fasta.gz'.format(url) for url in assemblies_ids]

start = time.time()


def download_assembly(url, file_name):
    """
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
            print('Retrying {0} ...{1}'.format(file_name, tries))

    return response


# We can use a with statement to ensure threads are cleaned up promptly
failures = []
with concurrent.futures.ThreadPoolExecutor(max_workers=30) as executor:
    # Start the load operations and mark each future with its URL
    for res in executor.map(download_assembly, ftp_urls, assemblies_ids):
        failures.append

print('Elapsed Time: %ss' % (time.time() - start))
