#! /usr/bin/env python3
"""
rename_htsf_fastq.py
Given a file of sample identifiers and fastq filenames, create a new naming scheme to use for symlinks.
"""

from __future__ import print_function

import os
import sys
import re
from collections import defaultdict

samples = defaultdict(dict)
for line in sys.stdin:

	if line.startswith("#"):
		continue

	iid, orig_fq = line.strip().split()[:2]

	prefix, fname = os.path.split(orig_fq)
	#m = re.search(r"_([ACGTN]+(?:_[A-Za-z0-9]+)*(?:_L[0-9]+))*_R([12])(?:_001)*\.fastq\.gz$", fname)
	#print(m.group(0,1,2))
	#run_id, which_mate = m.group(1,2)
	#which_mate = int(which_mate) - 1

	m = re.search(r"(.+)_R([12])(?:_001)*\.fastq\.gz$", fname)
	run_id, which_mate = m.group(1,2)
	which_mate = int(which_mate) - 1

	if not run_id in samples[iid]:
		samples[iid][run_id] = [None,None]
	samples[iid][run_id][which_mate] = orig_fq

prefix = "Pf"
ii = 0
for iid, lanes in samples.items():
	for lane, files in lanes.items():
		nfiles = len(list(filter(lambda x: x is not None, files)))
		new_id = "{}_{}{:04d}".format(iid, prefix, ii)
		new_files = [ os.path.join(iid, "{}_{}.fastq.gz".format(new_id, x)) for x in [1,2] ]
		print(iid, files[0], new_files[0], sep = "\t")
		print(iid, files[1], new_files[1], sep = "\t")
		print(iid, new_id, sep = "\t", file = sys.stderr)
		ii += 1
