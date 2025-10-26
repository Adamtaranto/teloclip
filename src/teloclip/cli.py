"""
Main CLI entry point for teloclip with sub-commands.
"""

import click

from teloclip._version import __version__


@click.group(
    help='A tool for the recovery of unassembled telomeres from soft-clipped read alignments.',
    invoke_without_command=True,
)
@click.version_option(version=__version__, prog_name='teloclip')
@click.pass_context
def main(ctx):
    """
    A tool for the recovery of unassembled telomeres from soft-clipped read alignments.

    Use sub-commands to filter alignments, extract reads, or extend contigs.

    Parameters
    ----------
    ctx : click.Context
        Click context object for passing information between commands.
    """
    # Ensure that ctx.obj exists and is a dict (in case `cli()` is called by itself)
    ctx.ensure_object(dict)

    # Check if no subcommand was invoked and show help
    if ctx.invoked_subcommand is None:
        click.echo(ctx.get_help())
        ctx.exit(0)


def register_commands():
    """Register sub-commands. Import here to avoid circular imports."""
    try:
        from teloclip.commands.extend import extend
        from teloclip.commands.extract import extract_cmd
        from teloclip.commands.filter import filter_cmd

        main.add_command(filter_cmd)
        main.add_command(extract_cmd)
        main.add_command(extend)
    except ImportError as e:
        # Handle gracefully during development
        click.echo(f'Warning: Could not import commands: {e}', err=True)


if __name__ == '__main__':
    # Register commands before running CLI
    register_commands()
    main()
else:
    # Register commands when module is imported
    register_commands()
