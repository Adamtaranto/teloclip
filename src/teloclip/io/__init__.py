"""
Reading and writing alignment and sequence files.

SAM text parsing, indexed BAM/FASTA access, buffered sequence output and input
format detection. Note that ``teloclip.io`` does not shadow the standard
library ``io`` module: Python 3 resolves bare ``import io`` absolutely, so
modules in this package still get the stdlib one.
"""
