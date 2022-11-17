#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------
This script accepts a list of accession numbers from the NCBI
databases and searches for linked identifiers to create a TSV
file with linked identifiers between the NCBI\'s databases.

Code documentation
------------------
"""


import os
import re
import csv
import time
import argparse

from Bio import Entrez


# regex expressions to identify identifier type
database_patterns = {'biosample': 'SAM[E|D|N][A-Z]?[0-9]+',
                     'bioproject': 'PRJ[E|D|N][A-Z][0-9]+',
                     'sra': '[E|D|S]RR[0-9]{6,}',
                     'refseq': 'GCF_[0-9]{9}.[0-9]+',
                     'genbank': 'GCA_[0-9]{9}.[0-9]+'}


def determine_id_type(identifier):
    """Determine the origin database for an accession number.

    Parameters
    ----------
    identifier : str
        Accession number from one of NCBI's databases.

    Returns
    -------
    match : str
        The name of the database the identifier belongs to.
        Returns None if it is not possible to identify the
        identifier type.
    """
    match = None
    for db, pat in database_patterns.items():
        db_match = re.findall(pat, identifier)
        if db_match != []:
            return db

    return match


def get_esearch_record(identifier, database):
    """Run an Entrez search and return the parsed results.

    Parameters
    ----------
    identifier : str
        Accession number from one of NCBI's databases.
    database : str
        NCBI database to query.

    Returns
    -------
    record : iter
        Multilevel data structure of Python lists and
        dictionaries with the query results, including
        the primary IDs linked to the accession number.
    """
    handle = Entrez.esearch(db=database, term=identifier)
    record = Entrez.read(handle)

    return record


def get_esummary_record(identifier, database):
    """Retrieve document summaries and return parsed results.

    Parameters
    ----------
    identifier : str
        Primary ID for a NCBI record.
    database : str
        NCBI database to query.

    Returns
    -------
    esummary_record : iter
        Multilevel data structure of Python lists and
        dictionaries with summary data about the record.
    """
    esummary_handle = Entrez.esummary(db=database, id=identifier, report='full')
    esummary_record = Entrez.read(esummary_handle, validate=False)

    return esummary_record


def get_elink_record(identifier, fromdb, todb):
    """Retrieve primary IDs for links to NCBI databases.

    Parameters
    ----------
    identifier : str
        Primary ID for a NCBI record.
    fromdb : str
        Database the input identifier belongs to.
    todb : str
        Database to search for linked identifiers.

    Returns
    -------
    elink_record : iter
        Multilevel data structure of Python lists and
        dictionaries with the query results, including
        primary IDs from `todb` linked to the input
        identifier.
    """
    elink_handle = Entrez.elink(dbfrom=fromdb, db=todb, id=identifier)
    elink_record = Entrez.read(elink_handle)

    return elink_record


def get_elink_id(elink_record):
    """Extract linked identifers from parsed results from Entrez.elink.

    Parameters
    ----------
    elink_record : iter
        Multilevel data structure of Python lists and
        dictionaries with the query results, including
        primary IDs from `todb` linked to the input
        identifier.

    Returns
    -------
    elink_ids : list
        List of primary IDs for one of NCBI's databases.
    """
    elink_id = elink_record[0]['LinkSetDb']
    if len(elink_id) > 0:
        elink_ids = elink_id[0]['Link']
        elink_ids = [i['Id'] for i in elink_ids]
    else:
        elink_ids = ''

    return elink_ids


def fetch_sra_accessions(identifiers):
    """Retrieve SRA accession numbers.

    Parameters
    ----------
    identifiers : list
        List of primary IDs for the SRA records.

    Returns
    -------
    sra_accessions : list
        List with SRA accession numbers.
    sequencing_platforms : list
        List with the sequencing platform for each
        SRA record.
    """
    sra_accessions = []
    sequencing_platforms = []
    for i in identifiers:
        # Get SRA summary
        sra_record = get_esummary_record(i, 'sra')

        # get SRA identifier
        sra_accession = re.findall(database_patterns['sra'],
                                   sra_record[0]['Runs'])
        if len(sra_accession) > 0:
            sra_accessions.append(sra_accession[0])

        # get sequencing platform
        sequencing_platform = sra_record[0]['ExpXml'].split('</Platform>')[0].split('>')[-1]
        sequencing_platforms.append(sequencing_platform)

    return [sra_accessions, sequencing_platforms]


def fetch_assembly_accessions(identifiers):
    """Retrieve Assembly accession numbers.

    Parameters
    ----------
    identifiers : list
        List with primary IDs for Assembly records.

    Returns
    -------
    refseq_accessions : list
        List with RefSeq accession numbers.
    genbank_accessions : list
        List with GenBank accession numbers.
    biosample_ids : list
        List with primary IDs for BioSample records
        linked to the Assembly records.
    biosample_accessions : list
        List with accession numbers for BioSample records
        linked to the Assembly records.
    """
    refseq_accessions = []
    genbank_accessions = []
    biosample_ids = []
    biosample_accessions = []
    for i in identifiers:
        # Get Assembly Summary
        assembly_record = get_esummary_record(i, 'assembly')
        # get RefSeq identifier
        refseq_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('RefSeq', '')
        if refseq_accession != '':
            refseq_accessions.append(refseq_accession)
        # get GenBank identifier
        genbank_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0]['Synonym'].get('Genbank', '')
        if genbank_accession != '':
            genbank_accessions.append(genbank_accession)

        # get Biosample accession number
        biosample_id = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleId', '')
        if biosample_id != '':
            biosample_ids.append(biosample_id)
        biosample_accession = assembly_record['DocumentSummarySet']['DocumentSummary'][0].get('BioSampleAccn', '')
        if biosample_accession != '':
            biosample_accessions.append(biosample_accession)

    return [refseq_accessions,
            genbank_accessions,
            list(set(biosample_ids)),
            list(set(biosample_accessions))]


# implement ID batch processing?
# (Bio.Entrez might support requests to get info for several IDs)
def main(input_file, output_file, email):

    # read identifiers
    with open(input_file, 'r') as infile:
        identifiers = infile.read().splitlines()

    # define email to make requests
    Entrez.email = email

    database_identifiers = {}
    # detect identifier type
    for i in identifiers:
        try:
            match = determine_id_type(i)
            if match is not None:
                print('{0:<} : {1}'.format(i, match.upper()))
                if match in ['refseq', 'genbank']:
                    match_db = 'assembly'
                else:
                    match_db = match
            else:
                print('{0:<} : {1}'.format('Could not determine database type.'))
                # process next identifier
                continue

            # get record data for identifier
            record = get_esearch_record(i, match_db)
            record_ids = record['IdList']

            # one Assembly identifier should only match one BioSample record?
            if match_db == 'assembly':
                refseq_accessions, genbank_accessions,\
                biosample_ids, biosample_accessions = fetch_assembly_accessions(record_ids)

                if len(biosample_ids) > 0:
                    # link to and get SRA accessions
                    elink_record = get_elink_record(biosample_ids[0], 'biosample', 'sra')
                    sra_ids = get_elink_id(elink_record)
                    sra_accessions, sequencing_platforms = fetch_sra_accessions(sra_ids)

                database_identifiers[i] = [refseq_accessions[0],
                                           genbank_accessions[0],
                                           biosample_accessions[0],
                                           ','.join(sra_accessions),
                                           ','.join(sequencing_platforms)]

            elif match_db == 'biosample':
                # link and get Assembly accessions
                elink_record = get_elink_record(record_ids[0], 'biosample', 'assembly')
                assembly_ids = get_elink_id(elink_record)
                refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

                # link to and get SRA accessions
                elink_record = get_elink_record(record_ids[0], 'biosample', 'sra')
                sra_ids = get_elink_id(elink_record)
                sra_accessions, sequencing_platforms = fetch_sra_accessions(sra_ids)

                database_identifiers[i] = [','.join(refseq_accessions),
                                           ','.join(genbank_accessions),
                                           i,
                                           ','.join(sra_accessions),
                                           ','.join(sequencing_platforms)]

            elif match_db == 'sra':
                # Get SRA Summary
                esummary_record = get_esummary_record(record_ids[0], 'sra')

                # get Biosample accession number (possible for one SRA Run accession to match multiple BioSample ids?)
                biosample_accession = esummary_record[0]['ExpXml'].split('<Biosample>')[-1].split('</Biosample>')[0]
                biosample_record = get_esearch_record(biosample_accession, 'biosample')
                biosample_id = biosample_record['IdList'][0]

                # find links to Assembly database
                # possible to get multiple RefSeq and GenBank ids
                elink_record = get_elink_record(biosample_id, 'biosample', 'assembly')
                assembly_ids = get_elink_id(elink_record)
                refseq_accessions, genbank_accessions = fetch_assembly_accessions(assembly_ids)[0:2]

                database_identifiers[i] = [','.join(refseq_accessions),
                                           ','.join(genbank_accessions),
                                           biosample_accession,
                                           i]
        except Exception:
            print('Could not retrieve data for {0}'.format(i))
            time.sleep(15)
            continue

        print('RefSeq: {0:<}\n'
              'GenBank: {1:<}\n'
              'BioSample: {2:<}\n'
              'SRA: {3:<}\n'
              'SequencingPlatform: {4:<}\n'.format(*database_identifiers[i]))

    # write output table
    output_header = 'InputID\tRefSeq\tGenBank\tBioSample\tSRA\tSequencingPlatform'
    output_lines = [output_header]
    output_lines.extend(['\t'.join([k]+v) for k, v in database_identifiers.items()])
    output_text = '\n'.join(output_lines)
    with open(output_file, 'w') as outfile:
        outfile.write(output_text+'\n')

    # delete .dtd files
    dtd_files = [file
                 for file in os.listdir(os.getcwd())
                 if file.endswith('.dtd')]
    for file in dtd_files:
        os.remove(file)


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', '--input-file', type=str,
                        required=True, dest='input_file',
                        help='Path to text file with NCBI database '
                             'identifiers (supports identifiers from '
                             'the Assembly, BioSample and SRA '
                             'databases).')

    parser.add_argument('-o', '--output-file', type=str,
                        required=True, dest='output_file',
                        help='Path to output TSV file with linked '
                             'identifiers between supported databases.')

    parser.add_argument('--email', type=str,
                        required=True, dest='email',
                        help='Email to perform requests with.')

    args = parser.parse_args()

    return args


if __name__ == '__main__':

    args = parse_arguments()
    main(**vars(args))
