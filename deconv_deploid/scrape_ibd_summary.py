#! /usr/bin/env python
"""
scrape_ibd_summary.py
Extract IBD state summaries from a dEploid log file.
"""

import os
import sys
import re
import argparse as ap

parser = ap.ArgumentParser(description = "Scrape IBD state summaries from dEploid log file.")
parser.add_argument("-i","--infile", type = ap.FileType("rU"),
			default = sys.stdin,
			help = "a dEploid log file" )
parser.add_argument("-s","--sample", type = str,
			default = None,
			help = "add this sample name as first column" )
args = parser.parse_args()

found_ibd = False

for line in args.infile:
	
	if re.match(r"\s*IBD probabilities", line):

		for ll in sys.stdin:

			this_line = ll.strip()
			if this_line == "":
				break

			combo, prop = this_line.split(":")
			prop = float(prop.strip())
			if args.sample is not None:
				print(args.sample, combo, prop, sep = "\t")
			else:
				print(combo, prop, sep = "\t")
			found_ibd = True

if not found_ibd:
	if args.sample is not None:
		print(args.sample, "0-0", 1.000, sep = "\t")
	else:
		print("0-0", 1.000, sep = "\t")
