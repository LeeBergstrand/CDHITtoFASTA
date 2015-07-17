#!/usr/bin/env python 
# Created by: Lee Bergstrand 
# Description: A provides functions for __main__.py
#
# ----------------------------------------------------------------------------------------
# ===========================================================================================================
# Imports:
import sys
import re
from Bio import SeqIO
from os import path

cluster_sequence_regex = re.compile("(\d*)\s*(\d*)aa,\s*(.*)\.\.\.\s*(?:at)?\s*([0-9\.\*]*)%?", re.IGNORECASE)


# ==========
# Classes:
# ==========
class Cluster(object):
	"""
	Definition of a CD-Hit cluster.

	:var cluster_id: The cluster's identity file.
	:var cluster_size: The sequence's accession.
	:var cluster_sequences: The cluster's sequences stored as a dictionary.
	"""

	def __init__(self, cluster_id, cluster_size, cluster_sequences=dict):
		self.cluster_id = int(cluster_id)
		self.cluster_sequences = dict(cluster_sequences)
		self.cluster_size = int(cluster_size)


# -------------------------------------------------------------------------------------------------
class ClusterSequence(object):
	"""
	Definition of a CD-Hit cluster sequence.

	:var seq_accession: The sequence's accession.
	:var seq_id: The sequence's identity within the cluster.
	:var seq_length: The length of the sequence.
	:var percent_identity: The sequences percent percent identity to the center sequence.
	:var center: Whether the sequence is the clusters center sequence.
	"""

	def __init__(self, seq_accession, seq_id, seq_length, percent_identity, center=False):
		self.seq_accession = str(seq_accession)
		self.seq_id = int(seq_id)
		self.seq_length = int(seq_length)
		self.percent_identity = float(percent_identity)
		self.center = center

	def is_center(self):
		if self.center:
			return True
		else:
			return False


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


# -----------------------------------------------------------------------------------------------------------
def parse_cd_hit_file(cd_cluster_string):
	# 2: Returns a list containing each cluster in the cd hit-file as a list of accession strings.
	cluster_object_list = []
	cluster_list = cd_cluster_string.split(">Cluster")
	del cluster_list[0]  # File starts with split character thus a empty element is created. Removes empty element.

	counter = 0
	for cluster in cluster_list:
		cluster_sequence_dict = {}

		cluster_rows = cluster.splitlines(True)
		del cluster_rows[0]

		for row in cluster_rows:
			cluster_sequence_object = parse_cluster_sequence(row)
			if cluster_sequence_object:
				cluster_sequence_dict[cluster_sequence_object.seq_accession] = cluster_sequence_object

		if cluster_sequence_dict:
			cluster_object_list.append(Cluster(counter, len(cluster_sequence_dict), cluster_sequence_dict))
		else:
			print("[Warning] Cluster " + counter + " is empty...")

		counter += 1

	return cluster_object_list


# -----------------------------------------------------------------------------------------------------------
def parse_cluster_sequence(cluster_string):
	"""
		Parses CD-Hit cluster sequence string into a cluster sequence object.

		:param cluster_string: The string representing a cluster sequence from a CD-Hit file.
		:return: A single cluster sequence object.
	"""
	sequence_features = cluster_sequence_regex.match(cluster_string)
	if sequence_features:
		seq_accession = sequence_features.group(3)
		seq_id = sequence_features.group(1)
		seq_length = sequence_features.group(2)
		percent_identity = sequence_features.group(4)
		center = False

		if percent_identity == '*':
			percent_identity = 100
			center = True

		cluster_sequence_object = ClusterSequence(seq_accession, seq_id, seq_length, percent_identity, center)
		return cluster_sequence_object
	else:
		print("[Warning] failed to parse: " + cluster_string)


# -----------------------------------------------------------------------------------------------------------
# 3: Extracts the accessions of each member of the CD-HIT cluster that contains a sequence from the sequence list file.
def getClusterAccessions(seqList, clusterList):
	clustersToCreateIntoFASTA = []
	for seq in seqList:
		for cluster in clusterList:
			if seq in cluster:
				if cluster not in clustersToCreateIntoFASTA:
					clustersToCreateIntoFASTA.append(cluster)
				break
	return clustersToCreateIntoFASTA


# -----------------------------------------------------------------------------------------------------------
# 4: Extracts the FASTA formatted sequece of each member of the CD-HIT cluster.
def getClusterFASTAs(clustersToCreateIntoFASTA, seqRecordDict):
	FASTAsToWrite = []
	for cluster in clustersToCreateIntoFASTA:
		ClusterFASTA = []
		for sequence in cluster:
			try:
				FASTA = seqRecordDict[sequence].format("fasta")
			except KeyError:
				print("Error: " + sequence + " could not be found in input FASTA file!")
				continue
			ClusterFASTA.append(FASTA)
		ClusterFASTAString = "".join(ClusterFASTA)
		FASTAsToWrite.append(ClusterFASTAString)
	return FASTAsToWrite


# ===========================================================================================================
"""# Main program code:
	
# House keeping...
argsCheck(4) # Checks if the number of arguments are correct.
AccessionExtractRegex = re.compile(">.*\.\.\.")

"""

"""

print ">> Extracting clusters which contain sequences from the sequence list."
ClusterAccessions = getClusterAccessions(seqList, clusterList)
print ">> Extracting clusters from the input FASTA file."
FASTAsToWrite = getClusterFASTAs(ClusterAccessions, seqRecords)

# Writes clusters to file.
print ">> Writing clusters to file in FASTA format."
ClusterCount = 1
for FASTA in  FASTAsToWrite:
	OutFile = seqFile.rstrip(".faa") + "Cluster" + str(ClusterCount) + ".faa"
	print ">> Writing Cluster " + str(ClusterCount) + " to file..."
	try:
		with open(OutFile,"w") as newFile:
			newFile.write(FASTA)
			newFile.close()
		print ">> Done."
	except IOError:
		print "Failed to write to " + OutFile
		sys.exit(1) # Aborts program. (exit(1) indicates that an error occurred)
	ClusterCount += 1
print ">> All Done."""
