#!/usr/bin/env python3
"""
Regenerate the coverage-decay figures embedded in the documentation.

Writes ``coverage-decay-analytic.png`` and ``coverage-decay-bootstrap.png``
into ``docs/images/`` (or a directory given with ``--out-dir``). The images
are committed artefacts, so run this after changing ``teloclip.model`` and
commit the result.

Requires the ``model`` extra: ``pip install -e '.[model]'``.

Usage::

    python scripts/generate_model_figures.py
    python scripts/generate_model_figures.py --out-dir /tmp/figs
"""

import argparse
import logging
from pathlib import Path

import matplotlib

# Agg before pyplot so the script runs headless (CI, ssh sessions).
matplotlib.use('Agg')

from teloclip.model.figures import generate_all  # noqa: E402

DEFAULT_OUT = Path(__file__).resolve().parent.parent / 'docs' / 'images'


def main(argv: 'list[str] | None' = None) -> None:
    """
    Parse arguments and regenerate the documentation figures.

    Parameters
    ----------
    argv : list of str or None, optional
        Argument list; None uses ``sys.argv``.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        '--out-dir',
        type=Path,
        default=DEFAULT_OUT,
        help=f'Directory to write the PNGs into (default: {DEFAULT_OUT})',
    )
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
    args.out_dir.mkdir(parents=True, exist_ok=True)
    for path in generate_all(str(args.out_dir)):
        logging.info('Wrote %s', path)


if __name__ == '__main__':
    main()
