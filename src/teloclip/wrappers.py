import re 
import os
from shlex import quote

def cleanID(s):
	"""Remove non alphanumeric characters from string. Replace whitespace with underscores."""
	s = re.sub(r"[^\w\s]", '', s)
	s = re.sub(r"\s+", '_', s)
	return s

def _minimap2_corrReads_cmd():
	#PB
	minimap2 -t num_threads -ax ava-pb --dual=yes reads reads > overlaps.sam 
	#ONT
	minimap2 -t num_threads -ax ava-ont --dual=yes reads reads > overlaps.sam 
	pass


def _minimap2_mapRef_cmd():
	minimap2 -ax map-pb ref.fa pacbio.fq.gz > aln.sam
	minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam
	pass



def _minimap2_prepMiniasm_cmd():
	minimap2 -x ava-pb -t8 pb-reads.fq pb-reads.fq | gzip -1 > reads.paf.gz


def _miniasm_cmd():
	
	return


def _racon_cor_cmd():
	racon -t num_threads -f reads overlaps.sam reads > polished_reads.fasta
	pass



######


def _nhmmer_command(exePath="nhmmer",modelPath=None,genome=None,evalue=None,nobias=False,matrix=None,cores=None,outdir=None):
	'''Construct the nhmmer command'''
	# Get model hmm basename
	model_base = os.path.splitext(os.path.basename(modelPath))[0]
	# Check for outdir
	if outdir:
		outdir = os.path.abspath(outdir)
		if not os.path.isdir(outdir):
				os.makedirs(outdir)
	else:
		outdir = os.getcwd()
	# Make subdir for nhmmer tab results
	outdir = os.path.join(outdir,"nhmmer_results")
	if not os.path.isdir(outdir):
		os.makedirs(outdir)
	# Create resultfile name
	outfile = os.path.join(os.path.abspath(outdir),model_base + ".tab")
	command = quote(exePath) + " --tblout " + quote(outfile)
	if cores:
			command += ' --cpu ' + str(cores)
	if evalue:
			command += ' -E ' + str(evalue)
	if nobias:
		command += ' --nobias'
	if matrix:
		command += ' --mxfile ' + quote(os.path.abspath(matrix))
	command += " --noali --notextw --dna --max " + quote(os.path.abspath(modelPath)) + " " + quote(os.path.abspath(genome))
	return command,outdir
