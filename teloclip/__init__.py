#!/usr/bin/env python
#python 3
#Contact, Adam Taranto, aptaranto@ucdavis.edu

######################################################################
# Filter SAM file for soft-clipped alignments containing unassembled #
# telomeric repeats.                                                 #
######################################################################

import os
import sys
import re

__version__ = "0.0.1"

class Error (Exception): pass

def isfile(path):
	path = os.path.abspath(path)
	if not os.path.isfile(path):
		print("Input file not found: %s" % path)
		sys.exit(1)
	else:
		return path

def read_fai(fai):
	path = isfile(fai)
	# Init empty dict
	ContigDict = dict()
	# Read fai_file to dict
	with open(path, "rU") as f:
		for line in f.readlines():
			li = line.strip().split()
			ContigDict[li[0]]=int(li[1])
	return ContigDict

def splitCIGAR(SAM_CIGAR):
	CIGARlist = list()
	for x in re.findall('[0-9]*[A-Z]',SAM_CIGAR):
		CIGARlist.append((re.findall('[0-9]*', x)[0],re.findall('[A-Z]', x)[0]))
	#174M76S --> [(174,M),(76,S)]
	#96S154M --> [(96,S),(154,M)]
	return CIGARlist

def checkClips(SAM_CIGAR):
	leftClipLen = None
	alnLen = None
	rightClipLen = None
	CIGARlist = splitCIGAR(SAM_CIGAR)
	# Check if first segment is soft-clipped
	if CIGARlist[0][1] == "S":
		leftClipLen = int(CIGARlist[0][0])
		if CIGARlist[1][1] == "M":
			alnLen= int(CIGARlist[1][0])
	# Check if last segment is soft-clipped
	if CIGARlist[-1][1] == "S":
		rightClipLen = int(CIGARlist[-1][0])
		if CIGARlist[-2][1] == "M":
			alnLen= int(CIGARlist[-2][0])
	return (leftClipLen, alnLen, rightClipLen)