#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
    
"""


import os
import argparse

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_str(string):
    """Reverse character order in input string.

    Parameters
    ----------
    string : str
        String to be reversed.

    Returns
    -------
    revstr : str
        Reverse of input string.
    """
    revstr = string[::-1]

    return revstr


def reverse_complement(dna_sequence):
    """Determine the reverse complement of given DNA strand.

    Parameters
    ----------
    strDNA : str
        String representing a DNA sequence.

    Returns
    -------
    revC_dna : str
        The reverse complement of the DNA sequence, without
        lowercase letters.

    Example:
        >>> reverse_complement('ATCGgcaNn')
        'NNTGCCGAT'
    """
    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    complement_bases = []
    for base in bases:
        if base in base_complement:
            complement_bases.append(base_complement[base])
        else:
            complement_bases.append(base.upper())

    complement_strand = ''.join(complement_bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def translate_sequence(dna_str, table_id):
    """Translate a DNA sequence using the BioPython package.

    Parameters
    ----------
    dna_str : str
        DNA sequence as string type.
    table_id : int
        Translation table identifier.

    Returns
    -------
    protseq : str
        Protein sequence created by translating the input
        DNA sequence.
    """
    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """Attempt to translate an input DNA sequence.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.
    method : str
        A string specifying the way the sequence will be
        oriented to attempt translation.
    table_id : int
        Translation table identifier.

    Returns
    -------
    List with following elements if translation is successful:
        protseq : str
            String representing the translated DNA sequence.
        myseq : str
            String representing the DNA sequence in the
            orientation used to translate it.
        Otherwise, returns string derived from captured exception.
    """
    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def check_str_alphabet(string, alphabet):
    """Determine if a string only has characters from specified alphabet.

    Parameters
    ----------
    string : str
        Input string.
    alphabet : str
        String that has all characters from desired alphabet.

    Returns
    -------
    "True" if sequence only has characters from specified
    alphabet and string "ambiguous or invalid characters" if
    it any of its characters is not in the alphabet.
    """
    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """Determine if length of input string is multiple of a specified number.

    Parameters
    ----------
    string : str
        Input string.
    number : int
        Integer that will be used to check if sequence length
        is multiple of.

    Returns
    -------
    "True" if the length of the sequence is a multiple of the
    specified number and "sequence length is not a multiple of number"
    if condition is not satisfied.
    """
    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id):
    """Check if sequence is valid and attempts to translate it.

    Calls several functions to ensure that the sequence only has
        'ACTG', is multiple of 3 and that it can be translated in any of 4
        different orientations. Stores exceptions so that it is possible to
        understand the sequence could not be translated.

    Parameters
    ----------
    dna_sequence : str
        String representing a DNA sequence.
    table_id (int):
        Translation table identifier.

    Returns
    -------
    If the sequence can be translated, a list with following elements:
        sequence : list
            A list with two elemets, the protein sequence and the DNA
            sequence in the correct orientation.
        coding_strand : str
            The strand orientation that had could be translated.
        Otherwise:
            exception_str : str
                A string containing the exceptions that determined
                that the sequence could not be translated.
    """
    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        sequence, exception_collector = retranslate(original_seq,
                                                    translating_methods[i],
                                                    table_id, strands[i],
                                                    exception_collector)

        i += 1
        if i == len(strands) or isinstance(sequence, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def retranslate(sequence, method, table_id, strands, exception_collector):
    """Attempt to translate a DNA sequence.

    Parameters
    ----------
    sequence : str
        String representing a DNA sequence.
    method : str
        A string specifying the sequence orientation that should be
        used to attempt translation.
    table_id : int
        Translation table identifier.
    strands : list
        List with 4 different orientations that can be checked.
    exception_collector : list
        List used to store all exceptions arising from translation
        attempts.

    Returns
    -------
    A list with following elements, if the sequence can be translated:
        translated_seq : list
            A list with the protein sequence and with the DNA sequence
            in the orientation used for translation.
        exception_collector : list
            A list with the exceptions that are captured when the
            sequence could not be translated.
        Otherwise:
            translated_seq : str
                A string with the exception/reason why the sequence
                could not be translated.
            exception_collector : list
                List with all exception that have been captured during
                translation attempts of the current sequence.
    """
    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def main(input_file, translation_table):

    parent_dir = os.path.dirname(input_file)
    basename = os.path.basename(input_file)
    file_prefix = basename.split('.fasta')[0]
    protein_file = os.path.join(parent_dir, file_prefix+'_protein.fasta')

    total = 0
    translated = 0
    with open(protein_file, 'a') as pf:
        for record in SeqIO.parse(input_file, 'fasta'):
            total += 1
            header = record.description
            dna = str(record.seq)
            translation_info = translate_dna(dna, translation_table)
            if isinstance(translation_info, str) is False:
                protein = str(translation_info[0][0])

                new_record = '>{0}\n{1}\n'.format(header, protein)

                pf.write(new_record)
                translated += 1
            else:
                print('{0}: {1}'.format(header, translation_info))

    if translated == 0:
        os.remove(protein_file)

    print('Translated {0}/{1}'.format(translated, total))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='Input FASTA file that contains DNA sequences '
                             'to be translated.')

    parser.add_argument('-t', '--translation-table', type=int,
                        required=False, default=11,
                        dest='translation_table',
                        help='Genetic code used to translate the DNA '
                             'sequences.')

    args = parser.parse_args()

    return args


if __name__ == "__main__":

    args = parse_arguments()
    main(**vars(args))
