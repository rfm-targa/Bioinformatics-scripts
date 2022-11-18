#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a BAM file and uses the SAMtools depth command
to create a file with the coverage depth for each nucleotide position.
The file with the coverage depth is sampled according to a user-defined
length interval and that information is processed to create a file with
the mean coverage depth for each extracted interval.

This script is written in Python 3 and does not require the installation
of additional Python packages (only requires SAMtools).

"""


import os
import argparse


def main(bam_file, output_directory, sample_interval):

    if not os.path.isfile(bam_file):
        raise Exception('\nInput file does not exist or cannot be '
                        'found in given path. Please provide a valid '
                        'path to a BAM file.')

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)
    else:
        print('\nOutput directory already exists. Moving on...')

    # run SAMtools depth command to compute coverage depth at each position
    print('\nUsing SAMtools to compute coverage depth '
          'at each nucleotide position...\n')
    depth_file = '{0}{1}'.format(os.path.basename(bam_file).split('.')[0],
                                 '.depth')
    depth_file_path = os.path.join(output_directory, depth_file)
    samtools_depth_cmd = 'samtools depth -aa {0} > {1}'.format(bam_file,
                                                               depth_file_path)
    os.system(samtools_depth_cmd)

    print('Importing coverage depth info and determining mean '
          'coverage depth for all intervals...\n')

    # read file with coverage depth for each nucleotide position
    with open(depth_file_path, 'r') as cov:
        bases_coverage = cov.readlines()
        bases_coverage = [line.strip().split('\t') for line in bases_coverage]

    # extract coverage information for consecutive intervals
    sample_info = []
    full_length_intervals = 0
    smaller_intervals = 0
    for c in range(0, len(bases_coverage), sample_interval):
        current_interval = bases_coverage[c:c+sample_interval]

        if len(current_interval) == sample_interval:
            full_length_intervals += 1
        else:
            smaller_intervals += 1

        # get identifer and position for starting sequence/contig
        start_position = current_interval[0][0:2]

        # get identifer and position for ending sequence/contig
        end_position = current_interval[-1][0:2]

        # calculate mean coverage value for current interval
        mean_coverage = sum([int(c[2])
                             for c in current_interval])/len(current_interval)

        # create line to add to report of mean coverage for each interval
        line = '{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(start_position[0],
                                                     start_position[1],
                                                     end_position[0],
                                                     end_position[1],
                                                     int(mean_coverage),
                                                     len(current_interval))

        sample_info.append(line)

    if sample_interval >= 1000:
        interval_str = str(sample_interval//1000)
    else:
        interval_str = str(sample_interval)

    output_file = '{0}{1}{2}{3}'.format(os.path.basename(bam_file.split('.')[0]), 
                                        '_', 
                                        interval_str, 
                                        'k_mean_depth.tsv')

    output_file_path = os.path.join(output_directory, output_file)

    output_header = ('Starting_sequence\tStarting_position\tEnding_sequence'
                     '\tEnding_position\tMean_coverage\tInterval_length\n')
    with open(output_file_path, 'w') as out:
        out.write(output_header)
        out.writelines('\n'.join(sample_info))

    print('Stored file with mean coverage length for {0} '
          'intervals of {1}bp and {2} smaller interval/s in '
          '{3}\n'.format(full_length_intervals, sample_interval,
                         smaller_intervals, output_file_path))

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str, required=True,
                        dest='input_file',
                        help='Path to the BAM file.')

    parser.add_argument('-o', '--output-directory', type=str, required=True,
                        dest='output_directory',
                        help='Path to the directory where files will be '
                             'stored.')

    parser.add_argument('-s', '--sample-interval', type=int, required=True,
                        dest='sample_interval',
                        help='The length interval value that will be used '
                             'to extract consecutive non-overlapping '
                             'regions to calculate the mean coverage.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
