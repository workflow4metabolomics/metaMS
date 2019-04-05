#!/usr/bin/env python
# vi: fdm=marker

import csv
import re
import argparse

# Get file cols {{{1
################################################################

def get_file_cols(file, preferred):

	cols = []

	with open(file if isinstance(file, str) else file.get_file_name(), 'r') as f:

		# Read file header
		reader = csv.reader(f, delimiter = "\t", quotechar='"')
		header = reader.next()

		preferred = preferred.split(',')

		# Determine default value
		perfect_matches = []
		partial_matches = []
		for p in preferred:
			for c in header:
				if c == p:
					perfect_matches.append(c) # Perfect match !
				elif re.match(p, c):
					partial_matches.append(c) # Keep this partial match in case we find no perfect match

		ordered_cols = perfect_matches + partial_matches
		for c in header:
			if not c in ordered_cols:
				ordered_cols.append(c)
		ordered_cols.append('NA')

		default = 0
		if len(perfect_matches) + len(partial_matches) == 0:
			default = len(ordered_cols) - 1

		# Build list of cols
		for i, c in enumerate(ordered_cols):
			cols.append( (c, c, i == default) )

	return cols

# Main {{{1
################################################################

if __name__ == '__main__':
    
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Script for getting column names in a csv file.')
    parser.add_argument('-f', help = 'CSV File (separator must be TAB)',       dest = 'file',    required = True)
    parser.add_argument('-p', help = 'List (comma separated values) of preferred column names for default one.',        dest = 'preferred',     required = True)
    args = parser.parse_args()
    args_dict = vars(args)
    
print(get_file_cols(**args_dict))