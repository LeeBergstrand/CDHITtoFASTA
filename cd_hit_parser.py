#!/usr/bin/env python

"""
Created by: Lee Bergstrand (Copyright 2015)

Description: A provides CD-Hit parsing functions and classes for lib.py and __main__.py
"""

# Imports & Setup:
import re
cluster_sequence_regex = re.compile("(\d*)\s*(\d*)aa,\s*>?(.*)\.\.\.\s*(?:at)?\s*([0-9\.\*]*)%?", re.IGNORECASE)


# ========
# Classes:
# ========

# -----------------------------------------------------------------------------------------------------------
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
def parse_cd_hit_file(cd_cluster_file_string):
	"""
	Parses the contents of a CD-Hit file into a list containing each CD-Hit cluster object.

	:param cd_cluster_file_string: String that contains the contents of a CD-Hit cluster file.
	:return: List of CD-Hit cluster objects.
	"""
	cluster_object_list = []
	cluster_list = cd_cluster_file_string.split(">Cluster")
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
