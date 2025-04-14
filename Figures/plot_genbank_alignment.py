#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script uses pyGenomeViz (https://github.com/moshi4/pyGenomeViz)
to create a figure representing the features and BLAST alignment
between a set of sequences in GenBank format.

For greater alignment sensitivity use the following BLASTn parameters:

'-word_size 28 -gapopen 10 -gapextend 5 -penalty -3 -reward 1 -xdrop_ungap 1 -xdrop_gap 1 -xdrop_gap_final 10'
"""

import os
import argparse

from pygenomeviz import GenomeViz
from pygenomeviz.parser import Genbank
from pygenomeviz.align import Blast, AlignCoord


def get_leftmost_alignments(alignment_coords):
	""" Determines the leftmost alignment for each query-reference pair.

	Parameters
	----------
	alignment_coords : list
		List of alignment coordinates (pyGenomeViz objects).

	Returns
	-------
	leftmost_values : dict
		A dictionary with the leftmost alignment start positions for each
		query-reference pair. The keys are the joined query and reference IDs,
		and the values are lists containing the start positions for the query
		and reference sequences.
	"""
	leftmost_values = {}
	# Extract start positions of the leftmost alignments
	for ac in alignment_coords:
		# Get query and reference IDs and start positions
		query_id = ac.query_id
		query_start = ac.query_start
		reference_id = ac.ref_id
		reference_start = ac.ref_start
		# Join query and reference IDs to use as key
		joined_id = '|'.join([query_id, reference_id])
		# If the joined ID is not in the dictionary, add it
		if joined_id not in leftmost_values:
			leftmost_values[joined_id] = [query_start, reference_start]
		# If the joined ID is in the dictionary, check if the current start position is smaller
		# than the stored one and update if necessary
		else:
			if query_start < leftmost_values[joined_id][0]:
				leftmost_values[joined_id] = [query_start, reference_start]

	return leftmost_values


def determine_offset(leftmost_coordinates):
	""" Determine the offset for each track based on the leftmost alignment for each sequence pair.

	Parameters
	----------
	"""
	# Store the offset for each track
	track_offsets = {}
	for k, v in leftmost_coordinates.items():
		# Split the joined ID to get the individual query and target IDs
		ids = k.split('|')
		# If the query leftmost position is greater than the target leftmost position
		if v[0] > v[1]:
			# Only determine offset for target if it still has none
			if ids[1] not in track_offsets:
				# Target offset is the difference between the query and target leftmost positions
				track_offsets[ids[1]] = v[0] - v[1]
				# Add zero offset for query if it is not in the dictionary
				if ids[0] not in track_offsets:
					track_offsets[ids[0]] = 0
				# If query already has a value, add its offset to the target offset
				else:
					track_offsets[ids[1]] += track_offsets[ids[0]]
		# If the target leftmost position is greater than the query leftmost position
		elif v[1] > v[0]:
			# Only determine offset for query if it still has none
			if ids[0] not in track_offsets:
				# Query offset is the difference between the target and query leftmost positions
				track_offsets[ids[0]] = v[1] - v[0]
				# Add zero offset for query if it is not in the dictionary
				if ids[1] not in track_offsets:
					track_offsets[ids[1]] = 0
				# If target already has a value, add its offset to the query offset
				else:
					track_offsets[ids[0]] += track_offsets[ids[1]]
		# If the leftmost positions are equal
		elif v[0] == v[1]:
			# If the query already has an offset value
			# Add the same offset value to the target
			if ids[0] in track_offsets:
				track_offsets[ids[1]] = track_offsets[ids[0]]
			# Set both to zero if neither has an offset
			else:
				track_offsets[ids[0]] = 0
				track_offsets[ids[1]] = 0

	return track_offsets


def main(input_files, output_directory, sequence_type, blast_options,
		 threads, minimum_length, minimum_identity, feature_color,
		 feature_style, label_type, alignment_color, alignment_inversed_color,
		 minimum_colorbar_value, track_alignment, colorbar_inverse, output_format):

	# List files in input directory
	gbk_list = list(map(Genbank, input_files))

	gv = GenomeViz()
	gv.set_scale_bar()

	# Run BLAST alignment before creating tracks to set offset if track_alignment == 'alignment'
	alignment_coords = Blast(gbk_list, seqtype=sequence_type, threads=threads, cmd_opts=blast_options).run()
	# Filter alignments by length and identity
	alignment_coords = AlignCoord.filter(alignment_coords, length_thr=minimum_length, identity_thr=minimum_identity)

	# Set track offset if user wants to align based on BLAST results
	# The 'alignment' option will align the tracks based on the leftmost BLAST alignment
	if track_alignment == 'alignment':
		leftmost_coordinates = get_leftmost_alignments(alignment_coords)
		# Determine offset for each track
		offsets = determine_offset(leftmost_coordinates)
	else:
		offsets = {gbk.name: track_alignment for gbk in gbk_list}

	# Create sequence tracks
	for gbk in gbk_list:
		track = gv.add_feature_track(gbk.name, gbk.get_seqid2size(), offset=offsets.get(gbk.name), align_label=False)
		# Add features to track
		for seqid, features in gbk.get_seqid2features("CDS").items():
			segment = track.get_segment(seqid)
			# lw = feature line width
			# fc = feature color
			# plotstyle = feature symbol
			segment.add_features(features, plotstyle=feature_style, label_type=label_type, fc=feature_color, lw=0.5)

	# Plot BLAST alignment links
	if len(alignment_coords) > 0:
		# Define minimum value for color scale
		# If minimum_colorbar_value is not None, use it
		# If minimum_colorbar_value is None, use the minimum identity of the alignments
		min_ident = minimum_colorbar_value if minimum_colorbar_value else int(min([ac.identity for ac in alignment_coords if ac.identity]))
		color, inverted_color = alignment_color, alignment_inversed_color
		for ac in alignment_coords:
			gv.add_link(ac.query_link, ac.ref_link, color=color, inverted_color=inverted_color, v=ac.identity, vmin=min_ident)
		# Add color bar
		if colorbar_inverse:
			gv.set_colorbar([color, inverted_color], vmin=min_ident)
		else:
			gv.set_colorbar([color], vmin=min_ident)

	if output_format == 'html':
		output_file = os.path.join(output_directory, "alignment_plot.html")
		gv.savefig_html(output_file)
	elif output_format == 'png':
		output_file = os.path.join(output_directory, "alignment_plot.png")
		gv.savefig(output_file, dpi=300)


def parse_arguments():

	parser = argparse.ArgumentParser(description="Generate alignment plots with pyGenomeViz.")

	parser.add_argument('-i', '--input-files', type=str, nargs='+',
						required=True, dest="input_files",
						help="Paths to the input files.")

	parser.add_argument('-o', '--output-directory', type=str,
						required=True, dest="output_directory",
						help="Path to the output directory.")

	parser.add_argument('-st', '--sequence-type', type=str,
						required=True, dest="sequence_type",
						choices=['nucleotide', 'protein'],
						help="Sequence type: nucleotide or protein.")

	parser.add_argument('-bo', '--blast-options', type=str,
						required=False, dest="blast_options",
						help="Additional BLASTn or tBLASTx options for alignment.")

	parser.add_argument('-t', '--threads', type=str,
						required=False, dest="threads",
						help="Number of threads passed to BLAST.")

	parser.add_argument('-ml', '--minimum-length', type=int,
						required=False, default=100,
						dest="minimum_length",
						help="Minimum length of an alignment to represent it.")

	parser.add_argument('-mi', '--minimum-identity', type=int,
						required=False, default=70,
						dest="minimum_identity",
						help="Minimum identity of an alignment to represent it.")

	parser.add_argument('-fc', '--feature-color', type=str,
						required=False, default='limegreen',
						dest="feature_color",
						help="Color of the features in the plot. "
							 "Any color in one of the formats accepted "
							 "by matplotlib can be used (https://matplotlib.org/stable/users/explain/colors/colors.html#color-formats).")

	parser.add_argument('-fs', '--feature-style', type=str,
						required=False, default='bigarrow',
						dest="feature_style",
						choices=['bigarrow', 'arrow', 'bigbox', 'box', 'bigrbox', 'rbox'],
						help="Style of the features in the plot.")

	parser.add_argument('-lt', '--label-type', type=str,
						required=False, dest="label_type",
						help="Attribute label value to use as feature label.")

	parser.add_argument('-ac', '--alignment-color', type=str,
						required=False, default='grey',
						dest="alignment_color",
						help="Color of the alignment links in the plot.")

	parser.add_argument('-aic', '--alignment-inversed-color', type=str,
						required=False, default='red',
						dest="alignment_inversed_color",
						help="Color of the alignment inversed links in the plot.")

	parser.add_argument('-mcv', '--minimum-colorbar-value', type=float,
						required=False,
						dest="minimum_colorbar_value",
						help="Minimum value for the colobar scale.")

	parser.add_argument('-ta', '--track-alignment', type=str,
					 	required=False, dest='track_alignment',
						choices=['left', 'center', 'right', 'alignment'],
						help='Define how to align tracks.')

	parser.add_argument('-ci', '--colorbar-inverse', action='store_true',
						required=False, default=False,
						dest="colorbar_inverse",
						help="If True, the color bar for inversions is displayed.")

	parser.add_argument('-of', '--output-format', type=str,
						required=False, default='html',
						dest="output_format",
						choices=['html', 'png'],
						help="Output file format. The defautl format is HTML. "
							 "The HTML allows to explore the plot interactively "
							 "and export in PNG and SVG formats.")

	return parser.parse_args()

if __name__ == "__main__":

	args = parse_arguments()
	main(**vars(args))
