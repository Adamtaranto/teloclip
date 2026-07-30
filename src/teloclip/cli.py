"""Main CLI entry point for teloclip with sub-commands."""

import importlib

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
    Teloclip is a tool for the recovery of unassembled telomeres from soft-clipped read alignments.

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
    """
    Register sub-commands.

    Each is imported independently so that one unavailable dependency does not
    take the others down with it. ``filter`` needs nothing beyond the standard
    library, click and rich, whereas ``extract`` and ``extend`` require pysam,
    pyfaidx and biopython. Importing all three together meant that in an
    environment without pysam no commands registered at all, including the one
    that would have worked.

    Imports are deferred to call time to avoid a circular import back to this
    module.
    """
    # (import path, attribute) for each sub-command.
    specs = [
        ('teloclip.commands.filter', 'filter_cmd'),
        ('teloclip.commands.extract', 'extract_cmd'),
        ('teloclip.commands.extend', 'extend'),
    ]

    for module_path, attribute in specs:
        try:
            module = importlib.import_module(module_path)
        except ImportError as e:
            click.echo(
                f'Warning: sub-command from {module_path} is unavailable: {e}',
                err=True,
            )
            continue

        main.add_command(getattr(module, attribute))


if __name__ == '__main__':
    # Register commands before running CLI
    register_commands()
    main()
else:
    # Register commands when module is imported
    register_commands()
