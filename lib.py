#!/usr/bin/env python

"""
Created by: Lee Bergstrand (Copyright 2015)
Description: A provides functions for __main__.py

Requirements:   - This software requires the Biopython module: http://biopython.org/wiki/Download
				- This software requires the cd_hit_parser.py module (included)
"""

# Imports:
import sys
from os import path

from Bio import SeqIO

from cd_hit_parser import *


# ==========
# Functions:
# ==========

# -----------------------------------------------------------------------------------------------------------
def check_file_extensions(reference_sequences_path, input_fasta_path, input_cluster_path):
	"""
	Checks file extensions for correctness.

	:param reference_sequences_path: The path for the reference sequences file.
	:param input_fasta_path: The path for the FASTA sequences file.
	:param input_cluster_path: The path for the CD-Hit cluster file.
	"""

	sequence_file_extension = path.splitext(reference_sequences_path)[-1]
	fasta_file_extension = path.splitext(input_fasta_path)[-1]
	cluster_file_extension = path.splitext(input_cluster_path)[-1]

	if not sequence_file_extension == ".txt":
		print("[Warning] " + sequence_file_extension + " may not be a txt file!")
	if not fasta_file_extension == ".faa":
		print("[Warning] " + fasta_file_extension + " may not be a FASTA file!")
	if not cluster_file_extension == ".clstr":
		print("[Warning] " + cluster_file_extension + " may not be a CD-Hit cluster file!")


# -----------------------------------------------------------------------------------------------------------
def get_reference_list(reference_sequences_path):
	"""
	Parses reference sequence file and returns a list of reference accessions.

	:param reference_sequences_path: The path for the reference sequences file.
	:return: List of reference sequence accession strings.
	"""

	try:
		infile = open(reference_sequences_path, "rU")
		sequences = infile.read()
		reference_sequences = sequences.splitlines()
		infile.close()
		return reference_sequences
	except IOError as e:
		print(str(e))
		sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# -----------------------------------------------------------------------------------------------------------
def get_fasta_dict(input_fasta_path):
	"""
	Reads a FASTA file and stores its sequences as a dictionary of Biopython sequence record objects.

	:param input_fasta_path: The path for the FASTA sequences file.
	:return: Dictionary of Biopython sequence record objects.
	"""

	try:
		new_file = open(input_fasta_path, "rU")
		sequence_record_dict = SeqIO.to_dict(SeqIO.parse(new_file, "fasta"))
		new_file.close()
		return sequence_record_dict
	except IOError as e:
		print(str(e))
		sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)


# -----------------------------------------------------------------------------------------------------------
def get_cluster_list(input_cluster_path):
	"""
	Reads a CD-HIT CLuster file passes its contents onto parse_cd_hit_file for parsing.

	:param input_cluster_path: The path for the CD-Hit cluster file.
	:return: List of Cluster objects.
	"""

	try:
		new_file = open(input_cluster_path, "rU")
		clustering = new_file.read()
		cluster_list = parse_cd_hit_file(clustering)
		new_file.close()
		return cluster_list
	except IOError as e:
		print(str(e))
		sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)
