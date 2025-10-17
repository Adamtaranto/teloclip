"""
Main CLI entry point for teloclip with sub-commands.
"""

import logging
import sys

import click

from teloclip._version import __version__
from teloclip.logs import init_logging


@click.group()
@click.option('--verbose', '-v', is_flag=True, help='Enable verbose logging')
@click.option('--quiet', '-q', is_flag=True, help='Suppress all but error messages')
@click.option(
    '--log-level',
    type=click.Choice(['DEBUG', 'INFO', 'WARNING', 'ERROR'], case_sensitive=False),
    help='Set specific log level',
)
@click.version_option(version=__version__, prog_name='teloclip')
@click.pass_context
def main(ctx, verbose, quiet, log_level):
    """
    A tool for the recovery of unassembled telomeres from soft-clipped read alignments.

    Use sub-commands to filter alignments, extract reads, or extend contigs.
    """
    # Ensure that ctx.obj exists and is a dict (in case `cli()` is called by itself)
    ctx.ensure_object(dict)

    # Configure logging based on options
    if log_level:
        level = getattr(logging, log_level.upper())
    elif verbose:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO

    # Initialize logging
    init_logging()
    logger = logging.getLogger()
    logger.setLevel(level)

    # Store logging level in context for sub-commands
    ctx.obj['log_level'] = level


def register_commands():
    """Register sub-commands. Import here to avoid circular imports."""
    try:
        from teloclip.commands.extract import extract_cmd
        from teloclip.commands.filter import filter_cmd

        main.add_command(filter_cmd)
        main.add_command(extract_cmd)
    except ImportError as e:
        # Handle gracefully during development
        click.echo(f'Warning: Could not import commands: {e}', err=True)


# Placeholder for extend command (will be implemented in Milestone 4)
@main.command('extend')
@click.pass_context
def extend_cmd(ctx):
    """Extend contigs using soft-clipped overhangs (Coming Soon)."""
    click.echo("The 'extend' command is not yet implemented.")
    click.echo('This feature will be available in a future release.')
    sys.exit(1)


if __name__ == '__main__':
    # Register commands before running CLI
    register_commands()
    main()
else:
    # Register commands when module is imported
    register_commands()
