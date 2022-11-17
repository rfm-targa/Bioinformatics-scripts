#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import os
import csv
import argparse


def main(input_files, output_directory, snps_table, ignore_terms):

    ignore_terms = [] if ignore_terms is None else ignore_terms    

    if os.path.isdir(output_directory) is False:
        os.mkdir(output_directory)
        print('Created {0}.'.format(output_directory))

    # read list of SNPs
    with open(snps_table, 'r') as infile:
        snps_list = list(csv.reader(infile, delimiter='\t'))

    # list of loci with SNPs
    # allow multiple possibilities for same position
    # e.g.: SNP or MNP or Complex
    stored = 0
    skipped = 0
    snps_pos = {}
    for l in snps_list[1:]:
        if any(t in l[1:4] for t in ignore_terms) is False:
            if l[0] not in snps_pos:
                snps_pos[l[0]] = {}
            snps_pos[l[0]].setdefault(l[5], [l[1:4]]).append(l[6])
            stored += 1
        else:
            skipped += 1

    print('Stored {0} relevant SNPs/indels and skipped {1} based on list of '
          'terms to ignore.'.format(stored, skipped))

    # list files in results folder
    results_files = [os.path.join(input_files, f)
                     for f in os.listdir(input_files)]

    print('Extracting SNPs in stored coordinates from {0} '
          'inputs...'.format(len(results_files)), end='')
    snps = {}
    for f in results_files:
        with open(f, 'r') as infile:
            current_results = list(csv.reader(infile, delimiter='\t'))

        # select SNPs based on position
        selected_snps = [l for l in current_results[1:] if l[1] in snps_pos]
        current_basename = os.path.basename(f)
        current_id = current_basename.split('.')[0]
        snps[current_id] = selected_snps

    # determine SNPs in positions
    snps_data = {}
    for k, v in snps.items():
        snps_data[k] = {}
        if len(v) > 0:
            for l in v:
                snps_data[k][l[1]] = (l[3], l[4])

        # add entry for missing SNPs
        for p in snps_pos:
            if p not in snps_data[k]:
                snps_data[k][p] = '-'

    # determine SNP presence
    snps_presence = {}
    for k, v in snps_data.items():
        snps_binary = {}
        for p, s in v.items():
            if s != '-' and s[0] in snps_pos[p]:
                relevant_snps = snps_pos[p][s[0]][1:]
                snps_binary[p] = 1 if s[1] in relevant_snps else 0
            else:
                snps_binary[p] = 0

        snps_presence[k] = snps_binary

    print('done.')

    # headers with loci names/product
    headers = ['FILE'] + ['{0}_{1}'.format(v[list(v.keys())[0]][0][0], k)
                          for k, v in snps_pos.items()] + ['Count']

    # SNPs list file
    snps_lines = []
    for k, v in snps_data.items():
        identifier = k
        current_snps = [v[p] for p in snps_pos]
        current_snps = ['{0}/{1}'.format(s[0], s[1]) if len(s) > 1 else s for s in current_snps]
        snps_sum = str(sum([1 for s in current_snps if s != '-']))
        current_snps.append(snps_sum)
        snps_lines.append([k]+current_snps)

    snps_lines = [headers] + snps_lines
    snps_lines = ['\t'.join(l) for l in snps_lines]
    snps_text = '\n'.join(snps_lines)

    snps_file = os.path.join(output_directory, 'snps_list.tsv')
    with open(snps_file, 'w') as outfile:
        outfile.write(snps_text+'\n')

    print('Wrote TSV file with the list of extracted SNPs/indels per input '
          'in start coordinates present in the SNPs table to '
          '{0}'.format(snps_file))

    # SNPs presence file
    presence_lines = []
    for k, v in snps_presence.items():
        identifier = k
        current_snps = [str(v[p]) for p in snps_pos]
        snps_sum = str(sum([1 for s in current_snps if s == '1']))
        current_snps.append(snps_sum)
        presence_lines.append([k]+current_snps)

    presence_lines = [headers] + presence_lines
    presence_lines = ['\t'.join(l) for l in presence_lines]
    presence_text = '\n'.join(presence_lines)

    presence_file = os.path.join(output_directory, 'snps_presence.tsv')
    with open(presence_file, 'w') as outfile:
        outfile.write(presence_text+'\n')

    print('Wrote TSV file with the presence absence matrix for extracted '
          'SNPs that match SNPs listed in the input table to '
          '{0}'.format(snps_file))

def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-files', type=str,
                        required=True, dest='input_files',
                        help='Path to directory that contains the TSV files '
                             'with the Snippy results, one file per strain.')

    parser.add_argument('-o', '--output-directory', type=str,
                        required=True, dest='output_directory',
                        help='Path to the output directory (will be created '
                             'if it does not exist).')

    parser.add_argument('-s', '--snps-table', type=str,
                        required=True, dest='snps_table',
                        help='Path to the TSV file that contains the list '
                             'of SNPs/indels.')

    parser.add_argument('-it', '--ignore-terms', nargs='+', type=str,
                        required=False, dest='ignore_terms',
                        help='List of terms to ignore. Lines in the SNPs '
                             'table that match any of the terms in the '
                             'Locus_tag, Gene or Product columns will be '
                             'skipped.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
