#!/usr/bin/env python

import os
import sys
import re
import shutil
import tempfile 
import subprocess
from itertools import groupby
from datetime import datetime
from collections import Counter

__version__ = "0.0.3"

###########

def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Error (Exception): pass

def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def _write_script(cmds,script):
    '''Write commands into a bash script'''
    f = open(script, 'w+')
    for cmd in cmds:
        print(cmd, file=f)
    f.close()

def syscall(cmd, verbose=False):
    '''Manage error handling when making syscalls'''
    if verbose:
        print('Running command:', cmd, flush=True)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        print('The following command failed with exit code', error.returncode, file=sys.stderr)
        print(cmd, file=sys.stderr)
        print('\nThe output was:\n', file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        raise Error('Error running command:', cmd)
    if verbose:
        print(decode(output))

def run_cmd(cmds,verbose=False,keeptemp=False):
    '''Write and excute script'''
    tmpdir = tempfile.mkdtemp(prefix='tmp.', dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = 'run_jobs.sh'
    _write_script(cmds,script)
    syscall('bash ' + script, verbose=verbose)
    os.chdir(original_dir)
    if not keeptemp:
        shutil.rmtree(tmpdir)


def getTimestring():
    """Return int only string of current datetime with milliseconds."""
    (dt, micro) = datetime.utcnow().strftime('%Y%m%d%H%M%S.%f').split('.')
    dt = "%s%03d" % (dt, int(micro) / 1000)
    return dt

def loadSam(samfile=None,contigs=None,maxBreak=50,minClip=1):
    SAM_QNAME   = 0
    SAM_RNAME   = 2
    SAM_POS     = 3
    SAM_CIGAR   = 5
    SAM_SEQ     = 9
    # Init dict
    alnDict = dict()
    # Add contig names as keys
    for name in contigs.keys():
        alnDict[name] = dict()
        alnDict[name]["L"] = list()
        alnDict[name]["R"] = list()
    # Read sam from stdin
    for line in samfile:
        # Skip header rows
        if line[0][0] == "@":
                continue
        samline = line.split('\t')
        # Check that aln contains soft clipping
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            # Get L/R clip lengths
            leftClipLen,rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (leftClipLen >= (int(samline[SAM_POS]) + minClip)):
                    # Overhang is on contig left
                    alnEnd = int(samline[SAM_POS]) + alnLen
                    try: alnDict[samline[SAM_RNAME]]["L"].append((samline[SAM_POS],alnEnd,leftClipLen,samline[SAM_SEQ],samline[SAM_QNAME]))
                    except: log("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
            # Check for right overhang
            if rightClipLen:
                try: ContigLen = contigs[str(samline[SAM_RNAME])]
                except: log("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
                alnEnd = int(samline[SAM_POS]) + alnLen 
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= maxBreak) and (alnEnd + rightClipLen >= ContigLen +1 ):
                    alnDict[samline[SAM_RNAME]]["R"].append((samline[SAM_POS],alnEnd,rightClipLen,samline[SAM_SEQ],samline[SAM_QNAME]))
    return alnDict
    

def StreamingSamFilter(samfile=None,contigs=None,maxBreak=50,minClip=1):
    '''Rewrite loadSam() as generator.'''
    SAM_QNAME   = 0
    SAM_RNAME   = 2
    SAM_POS     = 3
    SAM_CIGAR   = 5
    SAM_SEQ     = 9
    # Read sam from stdin
    for line in samfile:
        # Skip header rows
        if line[0][0] == "@":
                continue
        samline = line.split('\t')
        # Check that aln contains soft clipping
        if "S" in samline[SAM_CIGAR] and not "H" in samline[SAM_CIGAR]:
            # Get L/R clip lengths
            leftClipLen,rightClipLen = checkClips(samline[SAM_CIGAR])
            alnLen = lenCIGAR(samline[SAM_CIGAR])
            # Check for left overhang
            if leftClipLen:
                if (int(samline[SAM_POS]) <= maxBreak) and (leftClipLen >= (int(samline[SAM_POS]) + minClip)):
                    # Overhang is on contig left
                    alnEnd = int(samline[SAM_POS]) + alnLen
                    try: yield((samline[SAM_POS],alnEnd,leftClipLen,samline[SAM_SEQ],samline[SAM_QNAME],samline[SAM_RNAME],"L"))
                    except: log("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
            # Check for right overhang
            if rightClipLen:
                try: ContigLen = contigs[str(samline[SAM_RNAME])]
                except: log("Reference sequence not found in FAI file: " + str(samline[SAM_RNAME]))
                alnEnd = int(samline[SAM_POS]) + alnLen 
                # Check if overhang is on contig right end
                if ((ContigLen - alnEnd) <= maxBreak) and (alnEnd + rightClipLen >= ContigLen +1 ):
                    yield((samline[SAM_POS],alnEnd,rightClipLen,samline[SAM_SEQ],samline[SAM_QNAME],samline[SAM_RNAME],"R"))


def makeMask(killIdx,listlen):
    # makeMask([0,9], 10) =  [0,1,1,1,1,1,1,1,1,0]
    mask = [1 for i in range(listlen)]
    for x in killIdx:
        mask[x] = 0
    return mask

def filterList(data, exclude):
    # filterList([1,2,3,4,5,6,7,8,9,10],[0,9]) = [2,3,4,5,6,7,8,9]
    mask = makeMask(exclude,len(data))
    return (d for d,s in zip(data, mask) if s) 

def revComp(seq):
    """ Rev comp DNA string."""
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}[B] for B in x][::-1])
    return revcompl(seq)
    
def writeClip(idx,zpad,gap,seq,maplen):
    # leftpad idx ID
    padIdx = str(idx).zfill(zpad) + ':'
    # If gap between aln end(R) or start(L) and contig end, left pad softclip with '-'
    padseq = '-' * gap + seq
    # Format length of ref covered by alingment 
    readlen = 'LEN=' + str(maplen).rjust(6) 
    log('\t'.join([padIdx,readlen,padseq]))

def fasta2dict(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    contigDict = dict()
    for header in faiter:
        # Drop the ">"
        # Split on whitespace and take first item as name
        header = header.__next__()[1:].strip()
        name = header.split()[0]
        # Join all sequence lines to one.
        seq = "".join(s.strip() for s in faiter.__next__())
        contigDict[name] = (header, seq)
    return contigDict


def writefasta(outfile,name,seq,length=80):
    outfile.write('>' + str(name) + '\n')
    while len(seq) > 0:
        outfile.write(seq[:length] + '\n')
        seq = seq[length:]

def manageTemp(record=None, tempPath=None, scrub=False):
    """Create single sequence fasta files or scrub temp files."""
    if scrub and tempPath:
        try:
            os.remove(tempPath)
        except OSError:
            pass
    else:
        with open(tempPath, "w") as f:
            name,seq = record
            writefasta(f,name,seq,length=80)

def dochecks(args):
    """Housekeeping tasks: Create output files/dirs and temp dirs as required."""
    # Make outDir if does not exist else set to current dir.
    if args.temp:
        if not os.path.isdir(args.temp):
            os.makedirs(args.temp)
        tempParent = args.temp
    else:
        tempParent = os.getcwd()
    # Make temp directory
    os.makedirs(os.path.join(tempParent,"temp_" + getTimestring()))
    # Return full path to temp directory
    return tempDir


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []

def isfile(path):
    """
    Test for existence of input file.
    """
    if not os.path.isfile(path):
        print("Input file not found: %s" % path)
        sys.exit(1)
    else:
        return os.path.abspath(path)

def read_fai(fai):
    """
    Import fasta index file. Return dict of sequence names and lengths.
    """
    path = isfile(fai)
    # Init empty dict
    ContigDict = dict()
    # Read fai_file to dict
    with open(path, "rU") as f:
        for line in f.readlines():
            li = line.strip().split()
            ContigDict[li[0]]=int(li[1])
    return ContigDict

def addRevComplement(motifList):
    """
    Take list of DNA motif strings and return unique set of strings and their reverse complements.
    """
    revcompl = lambda x: ''.join([{'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}[B] for B in x][::-1])
    setList = list()
    for motif in motifList:
        setList.append(motif)
        setList.append(revcompl(motif))
    return set(setList)

def splitCIGAR(SAM_CIGAR):
    """
    Split CIGAR string into list of tuples with format (len,operator)
    """
    CIGARlist = list()
    for x in re.findall('[0-9]*[A-Z|=]',SAM_CIGAR):
        CIGARlist.append((int(re.findall('[0-9]*', x)[0]),re.findall('[A-Z]|=', x)[0]))
    #174M76S --> [(174,M),(76,S)]
    #96S154M --> [(96,S),(154,M)]
    return CIGARlist

def lenCIGAR(SAM_CIGAR):
    """
    Calculate length of alignment in reference sequence as sum of 
    match, read-deletion, splice, mismatch, and read-match block values.
    Ignore read-insertions, padding, hard and soft clip blocks.
    """
    alnLen = 0
    CIGARlist = splitCIGAR(SAM_CIGAR)
    for x in CIGARlist: # i.e. = [(174,M),(76,S)]
        if x[1] in set(['D','M','N','X','=']):
            alnLen += x[0]
    #Ignore operators in set('P','H','S','I')
    return alnLen

def checkClips(SAM_CIGAR):
    """
    Get lengths of soft-clipped blocks from either end of an alignment given a CIGAR string.
    """
    leftClipLen = None
    rightClipLen = None
    CIGARlist = splitCIGAR(SAM_CIGAR)
    # Check if first segment is soft-clipped
    if CIGARlist[0][1] == "S" :
        leftClipLen = int(CIGARlist[0][0])
    # Check if last segment is soft-clipped
    if CIGARlist[-1][1] == "S" :
        rightClipLen = int(CIGARlist[-1][0])
    return (leftClipLen,rightClipLen)

def crunchHomopolymers(motifList):
    """
    Take as input a list of target motifs, collapse poly-nucleotide tracks, return list of collapsed motifs.
    """
    # List to catch all collapsed motifs.
    crunchList = list()
    # For each motif
    for motif in motifList:
        # Create empty list to catch not repeated bases.
        noReps = list()
        # Walk through original motif base-by-base.
        for base in motif:
            # If list of kept bases in empty, add first base.
            if not noReps:
                noReps.append(base)
            # If list exists and base is not same as last, store new base.
            elif base != noReps[-1]:
                noReps.append(base)
        # Convert list to string and store new motif
        crunchList.append(''.join(noReps))
    # Convert to set to remove duplicates and return
    return list(set(crunchList))

def isClipMotif(samline,motifList,leftClip,rightClip,leftClipLen,rightClipLen,noPoly):
    """
    Extract terminal soft-clipped blocks from read sequence and test for presence of any DNA motif in motifList.
    """
    clipSeq = list()
    SAM_SEQ = 9
    if noPoly:
        if leftClip:
            clipSeq.append(crunchHomopolymers([samline[SAM_SEQ][0:leftClipLen]])[0])
        if rightClip:
            clipSeq.append(crunchHomopolymers([samline[SAM_SEQ][-rightClipLen:]])[0])
    else:
        if leftClip:
            clipSeq.append(samline[SAM_SEQ][0:leftClipLen])
        if rightClip:
            clipSeq.append(samline[SAM_SEQ][-rightClipLen:])
    # True if either clipped end sequence contains at least one instance of any motif.
    return any(s in x for s in motifList for x in clipSeq)


"""
CIGAR Operators
--------------------
D   Deletion; the nucleotide is present in the reference but not in the read
H   Hard Clipping; the clipped nucleotides are not present in the read.
I   Insertion; the nucleotide is present in the read  but not in the rference.
M   Match; can be either an alignment match or mismatch. The nucleotide is present in the reference.
N   Skipped region; a region of nucleotides is not present in the read
P   Padding; padded area in the read and not in the reference
S   Soft Clipping;  the clipped nucleotides are present in the read
X   Read Mismatch; the nucleotide is present in the reference
=   Read Match; the nucleotide is present in the reference

The output order in the array is “MIDNSHP=X” followed by a field for the NM tag. If the NM tag is not present, this field will always be 0.

M   BAM_CMATCH  0
I   BAM_CINS    1
D   BAM_CDEL    2
N   BAM_CREF_SKIP   3
S   BAM_CSOFT_CLIP  4
H   BAM_CHARD_CLIP  5
P   BAM_CPAD    6
=   BAM_CEQUAL  7
X   BAM_CDIFF   8
B   BAM_CBACK   9
NM  NM tag  10  
"""
