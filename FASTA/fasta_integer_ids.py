#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
Reads lines from a FASTA file given as input and substitutes
the identifiers of the headers by a single integer, starting
at the user defined value (default=1) and incrementing that
value by a step value (default=1) to assign a different integer
to the next sequence.

"""


import argparse


def main(input_fasta, output_fasta, starting_int, step):

    # read lines from FASTA file
    with open(input_fasta, 'r') as infile:
        lines = infile.readlines()

    # change sequence headers
    integer_id = starting_int
    for l in range(len(lines)):
        if lines[l].startswith('>'):
            lines[l] = '>{0}\n'.format(int(integer_id))
            integer_id += step

    # write lines with changed headers to output FASTA
    with open(output_fasta, 'w') as outfile:
        outfile.writelines(lines)

        sequence_count = sum([1 for line in lines if '>' in line])
        print('Wrote {0} sequences to {1}'.format(sequence_count,
                                                  output_fasta))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_fasta',
                        help='Path to the input FASTA file.')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output FASTA file.')

    parser.add_argument('-si', '--starting-integer', type=int,
                        required=False, default=1, dest='starting_integer',
                        help='Integer value that will be used to alter'
                             'the header of the first sequence.')

    parser.add_argument('-s', '--step', type=float, required=False,
                        dest='step', default=1,
                        help='The integer value used to alter headers'
                             'will be incremented by this value for each '
                             'sequence header that is altered.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
