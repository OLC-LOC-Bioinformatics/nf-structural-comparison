#!/usr/bin/env python

"""
Separate the contigs that were succesfully placed (scaffolded) from those that
were not, after scaffolding using RagTag.

Dependencies: biopython=1.79
"""

__author__ = 'Liam Brown'
__copyright__ = 'Crown Copyright 2022'
__license__ = 'GPL3'
__email__ = 'liam.brown@inspection.gc.ca'
__status__ = 'Prototype'

import os
import re
import sys
import argparse
from Bio import SeqIO

#-------------------------------------------------------------------------------
# Custom exceptions
#-------------------------------------------------------------------------------

class NoPlacedSeqsError(Exception):
    """
    Raise this error when no contigs from the draft genome could be placed
    (scaffoled) using the reference genome.
    """
    pass

#-------------------------------------------------------------------------------
# parse_arguments()
#-------------------------------------------------------------------------------

def parse_arguments():
    """
    Parse command-line arguments.
    """

    parser = argparse.ArgumentParser(
        description = 'Separate the contigs that were succesfully placed '
        '(scaffolded) from those that were not, after scaffolding using '
        'RagTag.')

    parser.add_argument('-s', '--scaffolded_draft', type = str, required = True,
        help = "Path to FASTA file produced by RagTag, containing the placed "
        "and unplaced sequences from the draft genome.")
    parser.add_argument('-p', '--placed', type = str,
        help = "Path to output FASTA file for saving placed sequences. "
        "Note: The '_RagTag' suffix will be removed from the header lines. "
        "Default: 'placed_seqs.fa'.",
        default = 'placed_seqs.fa')
    parser.add_argument('-u', '--unplaced', type = str,
        help = "Path to output FASTA file for saving placed sequences. "
        "Note: This file will not be created if all sequences were placed "
        "(scaffolded). "
        "Default: 'unplaced_seqs.fa'.",
        default = 'unplaced_seqs.fa')

    # If no arguments provided:
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    return args

#-------------------------------------------------------------------------------
# Other functions
#-------------------------------------------------------------------------------

def extract_seqs(scaffolded_draft, placed_fasta, unplaced_fasta):
    """
    Extract the placed (scaffolded) and unplaced sequences of the draft genome
    from RagTag's output FASTA file. Placed sequences are identified by FASTA
    record deflines terminating in '_RagTag'.

    :param scaffolded_draft (str): Path to FASTA file produced by RagTag, 
    containing the placed and unplaced sequences from the draft genome.
    :return placed_seqs (list): List of SeqRecord objects of placed sequences.
    :return unplaced_seqs (list): List of SeqRecord objects of unplaced 
    sequences.
    """
    placed_seqs = []
    unplaced_seqs = []

    if os.path.isfile(scaffolded_draft):

        if not os.path.isfile(placed_fasta) and not os.path.isfile(unplaced_fasta):

            for record in SeqIO.parse(scaffolded_draft, 'fasta'):
                # If the record.id ends with '_RagTag':
                if re.match('^.*_RagTag$', record.id):
                    # print(f'Record "{record.id}" from {os.path.basename(scaffolded_draft)} was placed.')
                    record.id = record.name = record.description = record.id.split('_RagTag')[0]
                    placed_seqs.append(record)
                else:
                    # print(f'Record "{record.id}" from {os.path.basename(scaffolded_draft)} was NOT placed.')
                    unplaced_seqs.append(record)

        elif os.path.isfile(placed_fasta):
            raise OSError(f"""
            File '{placed_fasta}' already exists. Please remove this file or
            move it to a different location and then re-run this script.
            """)

        elif os.path.isfile(unplaced_fasta):
            raise OSError(f"""
            File '{unplaced_fasta}' already exists. Please remove this file or
            move it to a different location and then re-run this script.
            """)
  
    return placed_seqs, unplaced_seqs

def write_output_files(placed_seqs, unplaced_seqs, scaffolded_draft, 
                       placed_fasta, unplaced_fasta):
    """
    Write output FASTA files containing placed (scaffolded) and unplaced
    sequences from the draft genome.

    :param placed_seqs (list): List of SeqRecord objects of placed sequences.
    :param unplaced_seqs (list): List of SeqRecord objects of unplaced 
    sequences.
    """
    # Write placed (scaffoled) and unplaced sequences to FASTA files
    if len(placed_seqs) > 0:
        print(f"Writing {len(placed_seqs)} placed sequences to '{placed_fasta}'.")
        SeqIO.write(placed_seqs, placed_fasta, 'fasta')
    else:
        raise NoPlacedSeqsError(f"""
        No placed sequences were detected in {scaffolded_draft}. Exiting...
        """)

    if len(unplaced_seqs) > 0:
        print(f"Writing {len(unplaced_seqs)} unplaced sequences to '{unplaced_fasta}'.")
        SeqIO.write(unplaced_seqs, unplaced_fasta, 'fasta')
    else:
        print(f"""
        All sequences were placed (scaffolded). Output '{unplaced_fasta}' was
        not created.
        """)

#-------------------------------------------------------------------------------
# main()
#-------------------------------------------------------------------------------

def main(args):

    placed_seqs, unplaced_seqs = extract_seqs(
        scaffolded_draft = args.scaffolded_draft,
        placed_fasta = args.placed,
        unplaced_fasta = args.unplaced)

    write_output_files(placed_seqs, unplaced_seqs,
        scaffolded_draft = args.scaffolded_draft,
        placed_fasta = args.placed,
        unplaced_fasta = args.unplaced)

#-------------------------------------------------------------------------------

if __name__ == '__main__':
    args = parse_arguments()
    main(args)