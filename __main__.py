#!/usr/bin/env python

"""
Created by: Lee Bergstrand (Copyright 2015)

Description:	Extracts CD-Hit clusters which contain reference proteins and stores them in FASTA format.

Requirements:	- This software requires the Biopython module: http://biopython.org/wiki/Download
				- This software requires the lib.py module (included)
				- This software requires the cd_hit_parser.py module (included)

usage: __main__.py [-h] [-i CLUSTER] [-s FASTA] [-r LIST]
"""

# Imports:
import argparse

from lib import *


# ----------------------------------------------------------------------------------------
def main(args):
	reference_list_path = args.reference_list[0]
	sequence_file_path = args.sequence_file[0]
	cluster_file_path = args.cluster_file[0]

	print("\nGenerating CD-HIT cluster FASTA files...\n")

	print("Sequence list file:  " + reference_list_path)
	print("Sequence FASTA file: " + sequence_file_path)
	print("CD-Hit cluster file: " + cluster_file_path + "\n")

	check_file_extensions(reference_list_path, sequence_file_path, cluster_file_path)

	reference_list = get_reference_list(reference_list_path)
	fasta_dict = get_fasta_dict(sequence_file_path)
	cluster_list = get_cluster_list(cluster_file_path)

	reference_clusters = []
	for cluster in cluster_list:
		for reference_accession in reference_list:
			if reference_accession in cluster.cluster_sequences.keys():
				reference_clusters.append(cluster)
				break
			else:
				continue

	fasta_list = []
	for cluster in reference_clusters:
		print("Cluster " + str(cluster.cluster_id) + " was found to have reference sequences.")
		cluster_accessions = cluster.cluster_sequences.keys()
		cluster_fastas = map(lambda x: fasta_dict[x].format("fasta"), cluster_accessions)
		cluster_tuple = (cluster.cluster_id, "".join(cluster_fastas))
		fasta_list.append(cluster_tuple)

	print("")
	for cluster_id, fasta in fasta_list:
		try:
			outfile = "Cluster" + str(cluster_id) + ".faa"
			print("Writing " + outfile)
			with open(outfile, "w") as new_file:
				new_file.write(fasta)
				new_file.close()
		except IOError as e:
			print(str(e))
			sys.exit(1)  # Aborts program. (exit(1) indicates that an error occurred)

	print("\nAll reference clusters written!")


# ----------------------------------------------------------------------------------------
if __name__ == '__main__':
	descriptor = """
	Extracts CD-Hit clusters which contain reference proteins and stores them in FASTA format.
	"""

	parser = argparse.ArgumentParser(description=descriptor)

	parser.add_argument('-i', '--cluster_file', metavar='CLUSTER', nargs=1, help='''
	CD-Hit cluster file which provides clustering information.''')

	parser.add_argument('-s', '--sequence_file', metavar='FASTA', nargs=1, help='''
	FASTA file which provides sequences to be extracted.''')

	parser.add_argument('-r', '--reference_list', metavar='LIST', nargs=1, help='''
	File of sequence identifiers (one per line) who's CD-HIT clusters should turned into FASTA files.''')

	cli_args = parser.parse_args()

	# At minimum we require a query, query DB and subject DB to proceed.
	proceed = True

	if cli_args.reference_list is None:
		print("Error: Missing sequence list path...")
		proceed = False

	if cli_args.sequence_file is None:
		print("Error: Missing sequence FASTA file path...")
		proceed = False

	if cli_args.cluster_file is None:
		print("Error: Missing query BLAST database path...")
		proceed = False

	if proceed:
		main(cli_args)
	else:
		print("")
		parser.print_help()
		print("")
