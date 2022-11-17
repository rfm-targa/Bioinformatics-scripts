#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script serves to analyse genome assemblies and create a report
with relevant statistics.

"""


import os
import bz2
import gzip
import time
import zipfile
import argparse
import itertools
import statistics as stats
import regex as re
from multiprocessing import Pool, cpu_count
import pandas as pd

from Bio import SeqIO
from Bio.SeqUtils import GC

COPEN = {
    "gz": gzip.open,
    "bz2": bz2.open,
    "zip": zipfile.ZipFile
}
MAGIC_DICT = {
    b"\x1f\x8b\x08": "gz",
    b"\x42\x5a\x68": "bz2",
    b"\x50\x4b\x03\x04": "zip"
}

"""
dict: Dictionary containing the binary signatures for three
compression formats (gzip, bzip2 and zip).
"""


def guess_file_compression(file_path, magic_dict=None):
    """Determine the type of compression of an input file.

    This function determines the compression of a given file
    by checking for a binary signature at the beginning of the file.
    These signatures are stored in the :py:data:`MAGIC_DICT`
    dictionary. The supported compression formats are gzip, bzip2
    and zip. If none of the signatures in this dictionary are found
    at the beginning of the file, it returns ``None``.

    Parameters
    ----------
    file_path : str
        Path to input file.
    magic_dict : dict, optional
        Dictionary containing the signatures of the compression types.
        The key should be the binary signature and the value should
        be the compression format. If left ``None``, it falls back to
        :py:data:`MAGIC_DICT`.

    Returns
    -------
    file_type : str or None
        If a compression type is detected, returns a string with
        the format. If not, returns ``None``.
    """
    if not magic_dict:
        magic_dict = MAGIC_DICT

    # only get the maximum number of characters that is needed to
    # guess any of the compressed formats

    max_len = max(len(x) for x in magic_dict)

    with open(file_path, "rb") as f:
        file_start = f.read(max_len)

    for magic, file_type in magic_dict.items():
        if file_start.startswith(magic):
            return file_type

    return None


def is_fasta(filename):
    """Check if a file is a FASTA file.

    Parameters
    ----------
    filename : str
        Full path to a FASTA file.

    Returns
    -------
    True if FASTA file, False otherwise.
    """
    try:
        with open(filename, "r") as handle:
            fasta = SeqIO.parse(handle, "fasta")

            # returns True if FASTA file, False otherwise
            return any(fasta)

    except UnicodeDecodeError:
        return False


def is_fasta_gz(filename):
    """Check if a file is a FASTA GZ file.

    Parameters
    ----------
    filename : str
        Full path to a FASTA file.

    Returns
    -------
    True if FASTA file, False otherwise.
    """
    with gzip.open(filename, "rt") as handle:
        fasta_gz = SeqIO.parse(handle, "fasta")

        # returns True if FASTA file, False otherwise
        return any(fasta_gz)


def check_if_list_or_folder(folder_or_list):
    """Check if a path represents a file or a directory.

    Parameters
    ----------
    folder_or_list : str
        Full path to the file or directory.

    Returns
    -------
    list_files : str
        If the input path represents a file, simply return
        input path. If input path represents a directory,
        lists files in the directory and stores paths to FASTA
        files in a file, returning the path to that file.
        Raises Exception if it is not possible to determine if
        the input path represents a file or a directory.
    """
    # check if input argument is a file or a directory
    if os.path.isfile(folder_or_list):
        list_files = folder_or_list

    elif os.path.isdir(folder_or_list):

        fasta_files = []

        for genome in os.listdir(folder_or_list):
            genepath = os.path.join(folder_or_list, genome)

            # do not continue if genepath is a dir
            if os.path.isdir(genepath):
                continue
            # check if file is a FASTA file
            if is_fasta(genepath):
                fasta_files.append(os.path.abspath(genepath))

            # check if file is a FASTA file
            elif is_fasta_gz(genepath):
                fasta_files.append(os.path.abspath(genepath))
        # if there are FASTA files
        if fasta_files:
            # store full paths to FASTA files
            with open('listGenes.txt', 'w') as f:
                for genome in fasta_files:
                    f.write(genome+'\n')
        else:
            raise Exception('There were no FASTA files in the given '
                            'directory. Please provide a directory '
                            'with FASTA files or a file with the list '
                            'of full paths to the FASTA files.')
        list_files = 'listGenes.txt'
    else:
        raise Exception('Input argument is not a valid directory '
                        'or file with a list of paths. Please provide'
                        'a valid input, either a folder with FASTA '
                        'files or a file with the list of full paths '
                        'to FASTA files (one per line).')

    return list_files


def verify_cpu_usage(cpu_to_use):
    """Verify CPU core count value.

    Parameters
    ----------
    cpu_to_use : int
        The number of CPU cores to use.

    Returns
    -------
    cpu_to_use : int
        The number of CPU cores adjusted to avoid using all
        available cores.
    """
    total_cpu = cpu_count()
    # do not allow a value of cpuToUse greater than the number of
    # cores/threads

    if cpu_to_use > total_cpu:
        """
        Detects present cpu and adjusts value so that machine doesn't crash
        """

        print('Warning! You have provided a CPU core count value that '
              'exceeds the number of cores in your machine!')
        print('Setting a different value for the CPU core count...')
        # define a value that is safe according to the number of
        # available cores/threads
        if total_cpu > 2:
            cpu_to_use = total_cpu - 2
        elif total_cpu == 2:
            cpu_to_use = 1

        print('CPU core count value set to: ', cpu_to_use)

    elif cpu_to_use < total_cpu and cpu_to_use > total_cpu - 2:

        print('Warning! You have provided a CPU core count value '
              'that is close to the maximum core count of your machine ('
              + str(cpu_to_use) + '/' + str(total_cpu) + '). This may '
              'affect your system responsiveness.')

    return cpu_to_use


def track_job(job, update_interval=3):
    """Track multiprocessing jobs."""
    while job._number_left > 0:
        remaining_tasks = job._number_left * job._chunksize
        print('Tasks remaining = {0}'.format(remaining_tasks))
        time.sleep(update_interval)


def flatten_list(list_to_flatten):
    """Flatten one level of a nested list.

    Parameters
    ----------
    list_to_flatten : list
        Input list to flatten.

    Returns
    -------
    flattened list : list
        Flattened list.
    """
    return list(itertools.chain(*list_to_flatten))


def analyse_report(report, nr_contigs, min_bp, max_bp, min_gc, max_gc, missing_data):
    """Determine if samples are of high quality based on provided thresholds.

    Parameters
    ----------
    report : dict
        Contains the assembly and annotation reports.
    nr_contigs : int
        Maximum number of contigs.
    min_bp : int
        Minimum number of total base pairs.
    max_bp : int
        Maximum number of total base pairs.
    min_gc : float
        Minimum GC content.
    max_gc : float
        Maximum GC content.

    Returns
    -------
    results : dict
        Dictionary with sample identifiers as keys and
        a list with the issues found for each sample as
        values.
    """
    for i, record in enumerate(report):
        current_results = []

        if record['Total_assembly_length'] < min_bp:
            current_results.append('Low_BP')

        elif record['Total_assembly_length'] > max_bp:
            current_results.append('High_BP')

        if record['GC_content'] < min_gc:
            current_results.append('Low_GC')

        elif record['GC_content'] > max_gc:
            current_results.append('High_GC')

        if record['Number_of_contigs'] > nr_contigs:
            current_results.append('Nr_contigs')

        if record['Missing_Data'] > missing_data:
            current_results.append('Too_many_N')

        report[i]['Warnings'] = ','.join(current_results)

    return report


def calc_n50(contig_sizes):
    """Calculate the N50 of an assembly.

    Parameters
    ----------
    contig_sizes : list
        List that contains the size of the contigs from
        an assembly.

    Returns
    -------
    l : int
        Calculated N50.
    """
    # Sort the contig sizes in descending order
    contig_sizes.sort(reverse=True)

    # Calculate n50
    total_bp = sum(contig_sizes)
    n50_threshold = total_bp * 0.5
    for l in contig_sizes:
        total_bp -= l
        if total_bp <= n50_threshold:
            # return the value of the largest contig equal
            # to or below n50 threshold
            return l


def analyse_assembly(assembly):
    """Analyse an assembly file.

    Parameters
    ----------
    assemblies : str
        Path to the assembly file in FASTA format.

    Returns
    -------
    results : dict
        Dictionary with key: value pairs for
        the sample identifier, number of contigs,
        average contig size, N50, total assembly length,
        GC content and missing data.
    """
    assembly_file = assembly

    # Get the sample name from the file
    sample = os.path.basename(assembly_file).split(".")[0]

    # Guess the file compression
    ftype = guess_file_compression(assembly_file)

    # File is not compressed
    if ftype is None:
        # Get the records of the assembly file
        records = list(SeqIO.parse(assembly_file, "fasta"))
    # This can guess the compression of gz, bz2 and zip.
    else:
        records = []
        with COPEN[ftype](assembly_file, "rt") as af:
            for record in SeqIO.parse(af, "fasta"):
                records.append(record)

    # Calculate the GC content
    all_gc_content = [GC(seq.seq) for seq in records]
    gc_content = stats.mean(all_gc_content) / 100

    # Get the total number of contigs in the assembly file
    nr_contigs = len(records)

    # Get the contig sizes
    sizes = [len(seq) for seq in records]

    # Calculate the average contig size
    avg_size = stats.mean(sizes)

    # Calculate the total assembly length
    total_length = sum(sizes)

    # Calculate the N50
    n50 = calc_n50(sizes)

    # Determine missing data
    missing_data = sum([rec.seq.count('N') for rec in records])

    n_blocks = []

    # Find all N blocks
    for rec in records:

        if 'N' in rec.seq != 0:
            n_blocks.append(list(re.findall('N+', str(rec.seq))))

    # N stats
    if len(n_blocks) == 0:
        num_blocks = 0
        min_missing_data = 0
        max_missing_data = 0
        median_missing_data = 0

    else:
        n_blocks = list(itertools.chain(*n_blocks))

        # real contigs
        num_blocks = len(n_blocks)

        # min missing data
        min_missing_data = len(min(n_blocks))

        # max missing data
        max_missing_data = len(max(n_blocks))

        # median value
        median_missing_data = stats.median(n_blocks.match)

    # Save the results in a dictionary
    results = {'Sample': sample,
               'Number_of_contigs': nr_contigs,
               'Average_contig_size': round(avg_size, 2),
               'N50': n50,
               'Total_assembly_length': total_length,
               'GC_content': round(gc_content, 3),
               'Missing_Data': missing_data,
               'number_of_N_blocks': num_blocks,
               'min_N': min_missing_data,
               'max_N': max_missing_data,
               'median_missing_data': median_missing_data}

    return results


def main(output_path, assembly_path, cpu, nr_contigs,
         minimum_number_of_bases, maximum_number_of_bases,
         minimum_gc_content, maximum_gc_content, missing_data):

    # avoid using all available cores
    # can lead to system unresponsiveness
    cpu_to_apply = verify_cpu_usage(cpu)
    # check if output directory exists

    if not os.path.exists(output_path):
        os.mkdir(output_path)

    if assembly_path:
        assemblies_file = check_if_list_or_folder(assembly_path)
        listGenes = []
        with open(assemblies_file, 'r') as gf:
            for gene in gf:
                gene = gene.rstrip('\n')
                listGenes.append(gene)
        listGenes.sort()
        os.remove('listGenes.txt')
        # List to save the results of the multiprocessing

        assembly_analysis_results = []

        print('Calculating assembly statistics...\n')

        p = Pool(processes=cpu_to_apply)
        r = p.map_async(analyse_assembly, listGenes,
                        callback=assembly_analysis_results.append)
        track_job(r)
        r.wait()

        print('\nAnalysing results...\n')

        # Flatten result nested list by one level
        results = flatten_list(assembly_analysis_results)

        # Analyse results
        results = analyse_report(results, nr_contigs,
                                 minimum_number_of_bases,
                                 maximum_number_of_bases,
                                 minimum_gc_content,
                                 maximum_gc_content,
                                 missing_data)

        # Print amount of Fails to console
        failed = sum([1 for k in results if len(k['Warnings']) > 0])
        print('The analysis detected {0} FAILS on a total '
              'of {1} assemblies. Check the report for more '
              'details.\n'.format(failed, len(listGenes)))
        print('Writing report...\n')

        # Convert dictionary into pandas DataFrame
        report = pd.DataFrame(results)

        # Write the final report
        output_report = os.path.join(output_path, 'final_report.tsv')
        report.to_csv(output_report, sep='\t',
                      encoding='utf-8', index=False)

    print('Execution Finished')


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-o', '--output', type=str, required=True,
                        dest='output_path', default=False,
                        help='Path to the output directory.')

    parser.add_argument('-i', '--assemblies', type=str, required=False,
                        dest='assembly_path', default=False,
                        help='Path to the directory containing the '
                             'geneome assemblies.')

    parser.add_argument('--cpu', type=int, required=False,
                        dest='cpu', default=2,
                        help='Number of CPU cores to use. Adjusts '
                             'the provided value if the it exceeds '
                             'the maximum number of available cores '
                             '(subtracts 2).')

    parser.add_argument('--nr_contigs', type=int, required=False,
                        dest='nr_contigs', default=350,
                        help='Maximum number of contigs allowed for '
                             'each assembly.')

    parser.add_argument('--min_bps', type=int, required=False,
                        dest='minimum_number_of_bases', default=1,
                        help='Minimum number of total bases '
                             'accepted for a genome assembly.')

    parser.add_argument('--max_bps', type=int, required=False,
                        dest='maximum_number_of_bases',
                        default=9999999999999999,
                        help='Maximum number of total bases '
                             'accepted for a genome assembly.')

    parser.add_argument('--min_gc', type=float, required=False,
                        dest='minimum_gc_content', default=0.0,
                        help='Minimum GC content value.')

    parser.add_argument('--max_gc', type=float, required=False,
                        dest='maximum_gc_content', default=1.0,
                        help='Maximum GC content value.')

    parser.add_argument('--min_N', type=int, required=False,
                        dest='missing_data', default=500,
                        help='Min number of N.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
