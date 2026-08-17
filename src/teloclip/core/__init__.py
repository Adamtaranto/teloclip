"""
Domain logic for teloclip: sequences, motifs, overhangs and extension.

Modules here describe *what* teloclip knows about terminal soft-clipped
alignments, independent of where the data came from or how it is presented.
They depend on ``teloclip.io`` only for streaming orchestration, and never on
``teloclip.report`` or ``teloclip.commands``.
"""
