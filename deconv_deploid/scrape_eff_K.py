#! /usr/bin/env python
"""
scrape_eff_K.py
Extract 'effective K' estimate from a dEploid log file.
"""

import os
import sys
import re
import argparse as ap

parser = ap.ArgumentParser(description = "Extract 'effective K' (weighted MOI) estimate from a dEploid log file.")
parser.add_argument("-i","--infile", type = ap.FileType("rU"),
			default = sys.stdin,
			help = "a dEploid log file" )
parser.add_argument("-s","--sample", type = str,
			default = None,
			help = "add this sample name as first column" )
args = parser.parse_args()

found_K = False

for line in args.infile:
	
	if re.match(r"\s*Effective_K", line):

		this_line = line.strip()
		_, eff_K = this_line.split(":")
		eff_K = float(eff_K.strip())

		next_line = next(args.infile).strip()
		_, hard_K = next_line.split(":")
		hard_K = int(hard_K.strip())

		if args.sample is not None:
			print(args.sample, hard_K, eff_K, sep = "\t")
		else:
			print(hard_K, eff_K, sep = "\t")
		found_K = True

if not found_K:
	if args.sample is not None:
		print(args.sample, -1, -1, sep = "\t")
	else:
		print(-1, -1, sep = "\t")
