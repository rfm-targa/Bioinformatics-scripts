#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION
    
    Accepts a matrix with results from the AlleleCall process of chewBBACA and 
    can apply two processes to that matrix: 
        
        - mask certain matrix elements by substituting those elements with other 
        characters; 
        - determine the number of shared loci and the number of allele differences 
        between each pair of samples represented in the AlleleCall matrix (with the 
        option to mask matrix elements before determining the number of shared alleles 
        and differences).
    
    The default masking option will substitute all ASM, ALM, NIPH, NIPHEM, PLOT3, 
    PLOT5, LOTSC and LNF cases with '0' and the 'INF-' prefix of inferred alleles 
    will always be removed to homogenize valid allele identifiers. 
    Passing a single word will change the default substitution value for that word. 
    To change specific matrix elements, the string should be formatted as:
        
                        'ASM=short,ALM=large,LNF=pop'
    
    which will change the default substitution value for ASM and ALM cases, maintaining 
    the default value of '0' to substitute remaining cases. The 'pop' expression serves 
    to signal that a type of case should not be substituted. 
    
    The ',' character is used to separate several substitution assignments and should only be 
    used that way.
    The '=' character is used to change the default substitution value for the right side 
    matrix element (ASM=short --> ASM cases will be substituted with the 'short' word) and 
    should only be used for that purpose.
    
    Execution examples:
        
        - 'mask' process:
            
            >>> python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> \\
                -p mask -o <masked_output_matrix> -mc ASM=short,ALM=large,LNF=pop
            
            Will substitute ASM elements with 'short' and ALM elements with 'large'.
            LNF cases will not be substituted due to the use of the 'pop' expression.
            NIPH, NIPHEM, PLOT3, PLOT5 and LOTSC cases will be substituted by the default
            '0' character. All 'INF-' prefixes will be removed from inferred alleles.
        
            >>> python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> \\
                -p mask -o <masked_output_matrix> -mc sub

            Will change the default subtitution value from '0' to 'sub' and all 
            ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5, LOTSC and LNF cases will be 
            substituted by 'sub'.
        
        - 'sims' process:
            
            >>> python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> \\
                -p sims -o <masked_output_matrix> -m True -mc ASM=short,LNF=pop
            
            Will substitute ASM elements with 'short'. LNF cases will not be substituted 
            due to the use of the 'pop' expression. ALM, NIPH, NIPHEM, PLOT3, PLOT5 and LOTSC 
            cases will be substituted by the default '0' character. All 'INF-' prefixes 
            will be removed from inferred alleles. After masking matrix values, the script 
            will determine a new matrix that has the number of loci shared between each pair 
            of samples above the diagonal and the number of allele differences below the diagonal.
            Masked values will not be considered when determining similarity, only cases where 
            each sample has a valid allele identifier will be considered.
            
            >>> python allelecall_simdiff_matrix.py -i <allele_call_matrix_file> \\
                -p sims -o <masked_output_matrix> -m False -mc ASM=short,LNF=pop
            
            Will not mask any matrix values due to the 'False' value passed to the 
            'm' argument. The number of sample similarities and differences will 
            be determined based on the full set of matrix elements because no matrix 
            value will be considered as masked.
"""


import csv
import argparse
import itertools


def mask_matrix(matrix_rows, masking_chars_dict):
    """ Masks matrix values. Masked values are substituted by either a default
        character or user-defined characters.
        
        Args:
            matrix_rows (list): list of sublists where each row corresponds has 
            the elements of one row of the original AlleleCall matrix.
            masking_chars_dict (dict): dictionary with matrix elements that should 
            be substituted as keys and the characters that will substitute those 
            elements as values.
            
        Returns:
            masked_matrix (list): list of sublists where each sublist has the elements 
            of one row of the matrix, masked according to the characters substitutions 
            specified in the 'masking_chars_dict' dictionary.
    """
    
    missing_data = masking_chars_dict
    
    # add first matrix row as header
    masked_matrix = []
    masked_matrix.append(matrix_rows[0])
    
    # for each matrix row that has an allelic profile
    for row in matrix_rows[1:]:
        # remove the 'INF-' prefix from each inferred allele to homogenize 
        # identifiers of valid alleles
        masked_row = [e.split('-')[1] if 'INF-' in e else e for e in row]
        
        # substitute the matrix elements that are keys in the missing_data dictionary
        masked_row = [missing_data[e] if e in missing_data else e for e in masked_row]
        
        # append each masked row
        masked_matrix.append(masked_row)

    return masked_matrix


def determine_simdiff(genomes_pairs, genomes_ids, masking_chars_dict):
    """ Determine the number of shared loci and the number of allele differences 
        between each possibe pair of samples.
        
        Args:
            genomes_pairs (list): list with pairs of genomes.
            genomes_ids (list): list with genomes identifiers.
            masking_chars_dict (dict): dictionary with matrix elements to 
            substitute as keys and the characters that should substitute those 
            elements as values.
                
        Returns:
            A list with two dictionaries:
                similarities (dict): dictionary with nested dictionaries. 
                Each nested dictionary stores the number of similarities between 
                a sample and all other samples.
                differences (dict): dictionary with nested dictionaries. 
                Each nested dictionary stores the number of differences between 
                a sample and all other samples.
    """

    # initialize empty dictionaries to store values
    similarities = {k:{} for k in genomes_ids}
    differences = {k:{} for k in genomes_ids}
    # iterate over genomes pairs
    invalid_values = list(masking_chars_dict.values())
    for pair in genomes_pairs:
        first_sample = pair[0]
        second_sample = pair[1]
        
        first_id = first_sample[0]
        second_id = second_sample[0]

        # determine valid comparisons
        pair_valid_comps =  [(i,j) for i, j in zip(first_sample[1:],second_sample[1:]) if i not in invalid_values and j not in invalid_values]
        # determine number of similarities
        pair_sims = sum([1 for t in pair_valid_comps if t[0] == t[1]])
        # subtract number of simiarities from total number of loci to get number of differences
        pair_diffs = len(pair_valid_comps) - pair_sims
        
        # store number of shared loci and number of allelic differences
        #similarities[first_id][second_id] = pair_sims
        similarities[first_id][second_id] = len(pair_valid_comps)
        differences[first_id][second_id] = pair_diffs

    return [similarities, differences]


def create_simdiff_matrix_rows(similarities, differences, genomes_ids):
    """ Creates a list with sublists where each sublist has the elements for one 
        row of a matrix with the number of shared loci between samples 
        above the diagonal and the number of allele differences below the diagonal.
        
        Args:
            similarities (dict): dictionary with nested dictionaries. 
            Each nested dictionary stores the number of similarities between 
            a sample and all other samples.
            differences (dict): dictionary with nested dictionaries. 
            Each nested dictionary stores the number of differences between 
            a sample and all other samples.
            genomes_ids (list): list with genomes identifiers.
            
        Returns:
            matrix_rows (list): list with sublists where each sublist has the 
            elements for one row of the matrix with the number of similarities 
            and differences between samples.
    """
    
    matrix_rows = []
    # add header with genomes identifiers as first line
    matrix_header = [''] + genomes_ids
    matrix_rows.append(matrix_header)
    # iterate over list of genomes identifiers
    for i, j in enumerate(genomes_ids):
        # get genome identifier index
        genome_index = i
        # get current genome identifier
        current_genome = j
        # get similarities and differences dictionaries
        current_sims = similarities[current_genome]
        current_diffs = differences[current_genome]
        
        # start appending line items by correct order
        # lines start with current_genome identifier
        new_row = []
        new_row.append(current_genome)
        # iterate over genomes identifiers list to add items to list in correct order
        for i2, j2 in enumerate(genomes_ids):
            
            # if the positions in the list are the same, it is the diagonal value
            if i2 == genome_index:
                new_row.append('0')
            # if the position is smaller than the position of the current genome
            # identifier in the list, we add values corresponding to differences
            # below the diagonal
            elif i2 < genome_index:
                diffs = str(current_diffs[j2]) if j2 in current_diffs else str(differences[j2][current_genome])
                new_row.append(diffs)
            # if the position is greater than the position of the current genome
            # identifier in the list, we add values corresponding to similarities
            # above the diagonal
            else:
                sims = str(current_sims[j2]) if j2 in current_sims else str(similarities[j2][current_genome])
                new_row.append(sims)

        matrix_rows.append(new_row)
    
    return matrix_rows


def sim_diff_matrix(matrix_rows, mask_missing, masking_chars_dict):
    """ Creates rows for a matrix that has the number of shared alleles between 
        pairs of samples above the diagonal and the number of differences below 
        the diagonal.
    
        Args:
            matrix_rows (list): list with sublists where each sublist has the 
            elements from one row of the original AlleleCall matrix.
            mask_missing (str): string indicating if the original AlleleCall 
            matrix values should be masked.
            masking_chars_dict (dict): dictionary that has the matrix elements 
            that should be masked as keys and the characters that should substitute 
            those matrix elements as values.
                
        Returns:
            matrix_rows (list): list of sublists where each sublist has the elements 
            of one row of the original matrix, masked or not masked depending on the 
            value of the 'mask_missing' argument.
            simdiff_matrix (list): list of sublists where each sublist has the number of 
            shared alleles between one sample and all other samples if the matrix values 
            are above the diagonal. Below the diagonal, values correspond to differences 
            between samples.
    """
    
    new_matrix_rows = matrix_rows
    
    if mask_missing == 'True':
        new_matrix_rows = mask_matrix(new_matrix_rows, masking_chars_dict)
    
    alleles = new_matrix_rows[1:]
    genomes_ids = [e[0] for e in alleles]

    # use itertools.combinations to find all pairs of samples
    genomes_pairs = list(itertools.combinations(alleles, 2))
    
    # iterate over pairs
    # find common elements by zipping lists and checking the number of common elements
    # in each tuple
    genomes_sims, genomes_diffs = determine_simdiff(genomes_pairs, genomes_ids, masking_chars_dict)

    # create lines for differences/similarities matrix
    simdiff_matrix = create_simdiff_matrix_rows(genomes_sims, genomes_diffs, genomes_ids)
    
    return [new_matrix_rows, simdiff_matrix]


def write_matrix(matrix_rows, output_matrix):
    """ Writes a matrix to a file.
    
        Args:
            matrix_lines (list): list of sublists where each sublist corresponds 
            to one row of the matrix.
            output_matrix (str): path to the file that should be created to 
            store the output matrix.
            
        Returns:
            Writes matrix rows (each sublist of the input list is a row) 
            to the output file.
    """
    
    # join matrix lines into chunk of text
    concat_lines = ['\t'.join(line) for line in matrix_rows]
    matrix_text = '\n'.join(concat_lines)
    
    # write matrix to output file
    with open(output_matrix, 'w') as output_file:
        output_file.writelines(matrix_text)


def main(allele_call_matrix, process_name, output_matrix, mask_missing, masking_character):
    """ Accepts a matrix with results from the AlleleCall process of chewBBACA and 
        can apply two processes to that matrix, mask certain matrix elements by 
        substituting those elements with other characters or determine the number 
        of shared alleles and differences between each pair of samples represented in 
        the AlleleCall matrix (with the option to mask matrix elements before determining 
        the number of shared alleles and differences).
        
        Args:
            allele_call_matrix (str): path to the file with the AlleleCall matrix.
            process_name (str): name of the process to be applied, either 'mask' or 'sims'.
            output_matrix (str): path to the file that will be created to store the output matrix.
            mask_missing (str): string to signal if matrix values should be masked or not, either
            'True' or 'False'. Only matters for the 'sims' process, in the 'mask' process matrix 
            elements will be masked independently of this argument value.
            masking_character (str): the character or characters that should be used to substitute 
            matrix elements. The default option will substitute all ASM, ALM, NIPH, NIPHEM, PLOT3, 
            PLOT5, LOTSC, LNF with '0' and the 'INF-' prefix of inferred alleles will always be removed 
            to homogenize valid alleles identifiers. Passing a single word will change the default 
            substitution value to that word. To change specific matrix elements, the string should be 
            formatted as 'ASM=short,ALM=large,LNF=pop', which will change the default substitution 
            value for ASM and ALM cases, maintaining the default value of '0' to substitute remaining 
            cases. The 'pop' expression serves to signal that a type of case should not be substituted.
            
        Returns:
            Outputs a file to the path specified in 'output_matrix' with a matrix that will correspond 
            to the AlleleCall matrix with masked elements if the user selected the 'mask' process or with 
            the number of shared alleles and differences if the user selected the 'sims' process.
    """

    # import AlleleCall matrix
    with open(allele_call_matrix, 'r') as matrix:
        matrix_lines = list(csv.reader(matrix, delimiter='\t'))

    # create substitutions dictionary
    masking_dict = {'ALM':'0', 'ASM':'0', 'LNF':'0', 'NIPH':'0', 'NIPHEM':'0', 
                    'PLOT3':'0', 'PLOT5':'0', 'LOTSC':'0'}

    # alter dictionary when users want to alter dictionary values
    # for a value different than default '0' or when cases to substitute
    # should not be all substituted by the same character
    chars = masking_character.split(',')
    if '=' in chars[0]:
        for char in chars:
            if '=' in char:
                c = char.split('=')
                if c[1] != 'pop':
                    masking_dict[c[0]] = c[1]
                else:
                    masking_dict.pop(c[0])
                
                try:
                    new_char = int(c[1])
                    if new_char != '0':
                        print('WARNING: You have chosen {0} to substitute for {1}. This might conflict '
                              'with valid allele identifiers.'.format(c[1], c[0]))
                except:
                    continue

    # alter dictionary values when users do not want to substitute all cases with '0'
    elif len(chars) == 1 and '=' not in chars[0]:

        masking_dict = {k:chars[0] for k, v in masking_dict.items()}

        try:
            new_char = int(chars[0])
            print('WARNING: You have chosen {0} to substitute all cases. This might conflict '
                  'with valid allele identifiers.'.format(new_char))
        except:
            pass

    # simply mask certain matrix values
    if process_name == 'mask':
        masked_matrix = mask_matrix(matrix_lines, masking_dict)
        write_matrix(masked_matrix, output_matrix)

    # Determine the number of pairwise similarities and differences
    # with the option to mask matrix values before
    elif process_name == 'sims':
        masked_matrix, simdiff_matrix = sim_diff_matrix(matrix_lines, mask_missing, masking_dict)
        write_matrix(masked_matrix, 'masked_'+output_matrix)
        write_matrix(simdiff_matrix, output_matrix)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--allele_call_matrix', type=str, required=True, dest='allele_call_matrix',
                        help='Path to the file with the AlleleCall matrix (default name given by chewBBACA is results_alleles.tsv).')

    parser.add_argument('-p', '--process_name', type=str, required=True, dest='process_name', choices=['mask', 'sims'],
                        help='Choose the type of process to run. The input matrix can be masked (mask mode) to substitute ASM, ALM, NIPH, NIPHEM, '
                             'PLOT3, PLOT5, LOTSC, LNF with given characters. The "INF-" prefix of inferred alleles is removed from the '
                             ' string by default to homogenize valid allele identifiers. '
                             'Another option that is available (mode sims) is the creation of a matrix with the number of shared alleles between pairs '
                             'of genomes/assemblies above the diagonal and with the number of differences below the diagonal (with the option '
                             'to mask matrix values before determining similarities and differences).')

    parser.add_argument('-o', '--output_matrix_filename', type=str, required=True, dest='output_matrix_filename',
                        help='Path to the output file that will be created to store the new matrix.')

    parser.add_argument('-m', '--mask_missing_data', type=str, required=False, dest='mask_missing_data', choices=['True', 'False'], default='True',
                        help='Define if the AlleleCall results corresponding to ASM, ALM, NIPH, NIPHEM, PLOT3, PLOT5, '
                             'LOTSC, LNF and the prefix "INF-" of inferred alleles should be masked (default=True). '
                             'This argument is only necessary for the "sims" process.')

    parser.add_argument('-mc', '--masking_character', type=str, required=False, dest='masking_character', default='0',
                        help='Define the character that will substitute the cases to be masked. By default, all cases '
                             'that are not valid alleles are substituted by "0" and inferred alleles are stripped from '
                             'the "INF-" prefix. Defining another single character will substitute all cases with that character. '
                             'Different cases can be substituted by different characters but the characters must be given in '
                             'the format "ASM=short,ALM=long,NIPH=para,NIPHEM=parexc". Cases that are not specified will be '
                             'substituted by the default and the value "pop" will remove the specified cases from the list of cases '
                             'that should be substituted.')

    args = parser.parse_args()

    return [args.allele_call_matrix, args.process_name, args.output_matrix_filename, 
            args.mask_missing_data, args.masking_character]


if __name__ == "__main__":

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3], args[4])

