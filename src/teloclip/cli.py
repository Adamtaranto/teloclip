"""Main CLI entry point for teloclip with sub-commands."""

import importlib

import click

from teloclip._version import __version__


class _LazyGroup(click.Group):
    """
    Click group that registers its sub-commands on first use.

    Registration used to run as an import side effect, which meant that merely
    importing :mod:`teloclip.cli` pulled in pysam, pyfaidx and biopython and
    mutated the group. Deferring it to the first ``get_command`` or
    ``list_commands`` keeps the import cheap and side-effect free while leaving
    the observable CLI behaviour unchanged.
    """

    # Class-level default. _ensure_registered shadows it with an instance
    # attribute, which avoids overriding __init__ just to set a flag.
    _registered = False

    def _ensure_registered(self) -> None:
        """Import and attach the sub-commands, at most once."""
        if not self._registered:
            # Set first: register_commands() adds to this group, and a failure
            # part way through must not cause a second attempt to re-add the
            # commands that did succeed.
            self._registered = True
            register_commands()

    def get_command(self, ctx, cmd_name):
        """
        Look up a sub-command by name, registering the set if needed.

        Parameters
        ----------
        ctx : click.Context
            Click context object.
        cmd_name : str
            Name of the sub-command being resolved.

        Returns
        -------
        click.Command or None
            The command, or None if no sub-command has that name.
        """
        self._ensure_registered()
        return super().get_command(ctx, cmd_name)

    def list_commands(self, ctx):
        """
        List sub-command names, registering the set if needed.

        Parameters
        ----------
        ctx : click.Context
            Click context object.

        Returns
        -------
        list of str
            Sorted sub-command names.
        """
        self._ensure_registered()
        return super().list_commands(ctx)


@click.group(
    cls=_LazyGroup,
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
    main()
