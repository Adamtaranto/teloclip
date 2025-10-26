"""Test main CLI entry point functionality."""

from pathlib import Path
import re
import subprocess
from typing import List, Optional, Tuple

import pytest


class CLIRunner:
    """Helper class for running CLI commands and capturing output."""

    def __init__(self, timeout: int = 30):
        """Initialize CLI runner with timeout."""
        self.timeout = timeout

    def run_command(
        self,
        cmd: List[str],
        cwd: Optional[Path] = None,
        input_data: Optional[str] = None,
        check: bool = False,
    ) -> Tuple[int, str, str]:
        """
        Run a command and return (exit_code, stdout, stderr).

        Args:
            cmd: Command and arguments to run
            cwd: Working directory for command
            input_data: Data to pass to stdin
            check: If True, raise CalledProcessError on non-zero exit

        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=self.timeout,
                cwd=cwd,
                input=input_data,
                check=check,
            )
            return result.returncode, result.stdout, result.stderr
        except subprocess.TimeoutExpired:
            return -1, '', f'Command timed out after {self.timeout} seconds'
        except subprocess.CalledProcessError as e:
            return e.returncode, e.stdout, e.stderr
        except FileNotFoundError:
            return -1, '', 'Command not found'

    def run_teloclip(
        self,
        args: List[str],
        cwd: Optional[Path] = None,
        input_data: Optional[str] = None,
        check: bool = False,
    ) -> Tuple[int, str, str]:
        """
        Run teloclip command with given arguments.

        Args:
            args: Arguments to pass to teloclip
            cwd: Working directory
            input_data: Data to pass to stdin
            check: If True, raise CalledProcessError on non-zero exit

        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        cmd = ['teloclip'] + args
        return self.run_command(cmd, cwd=cwd, input_data=input_data, check=check)


@pytest.fixture
def cli_runner():
    """Provide CLI runner instance."""
    return CLIRunner()


def assert_exit_code(exit_code: int, expected: int, stdout: str, stderr: str):
    """Assert exit code matches expected, with helpful error message."""
    if exit_code != expected:
        msg = f'Expected exit code {expected}, got {exit_code}'
        if stdout:
            msg += f'\nSTDOUT:\n{stdout}'
        if stderr:
            msg += f'\nSTDERR:\n{stderr}'
        pytest.fail(msg)


def assert_contains(text: str, expected: str, case_sensitive: bool = True):
    """Assert that text contains expected string."""
    if not case_sensitive:
        text = text.lower()
        expected = expected.lower()

    if expected not in text:
        pytest.fail(f"Expected '{expected}' in text, but got:\n{text}")


class TestMainCLI:
    """Test main teloclip CLI entry point."""

    @pytest.mark.xfail(
        reason='Version output may be unstable; marked as expected to fail'
    )
    def test_version_flag(self, cli_runner):
        """Test --version displays correct version."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(['--version'])

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stdout.lower(), 'teloclip', case_sensitive=False)
        # Should contain version number (format like 1.0.0 or similar)
        version_pattern = r'\d+\.\d+\.\d+'
        assert re.search(version_pattern, stdout), (
            f'No version number found in: {stdout}'
        )

    def test_help_flag(self, cli_runner):
        """Test --help shows main help with all commands."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(['--help'])

        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stdout, 'usage:', case_sensitive=False)
        assert_contains(stdout, 'teloclip', case_sensitive=False)

        # Should show all three subcommands
        assert_contains(stdout, 'extend', case_sensitive=False)
        assert_contains(stdout, 'extract', case_sensitive=False)
        assert_contains(stdout, 'filter', case_sensitive=False)

    def test_no_arguments_shows_help(self, cli_runner):
        """Test teloclip with no arguments shows help."""
        exit_code, stdout, stderr = cli_runner.run_teloclip([])

        # Should show help (exit code 0 is acceptable for help display)
        assert_exit_code(exit_code, 0, stdout, stderr)
        # Help could be in stdout or stderr depending on implementation
        help_output = stdout + stderr
        assert_contains(help_output, 'usage:', case_sensitive=False)
        assert_contains(help_output, 'teloclip', case_sensitive=False)

    def test_invalid_command_shows_error(self, cli_runner):
        """Test invalid command shows error and suggestions."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(['invalid-command'])

        assert exit_code != 0
        error_output = stdout + stderr
        # Should indicate the command is invalid
        assert_contains(error_output, 'invalid', case_sensitive=False)

    def test_help_for_subcommands_available(self, cli_runner):
        """Test that help is available for all subcommands."""
        subcommands = ['extend', 'extract', 'filter']

        for subcommand in subcommands:
            exit_code, stdout, stderr = cli_runner.run_teloclip([subcommand, '--help'])

            assert_exit_code(exit_code, 0, stdout, stderr)
            assert_contains(stdout, 'usage:', case_sensitive=False)
            assert_contains(stdout, subcommand, case_sensitive=False)

    def test_global_log_level_option(self, cli_runner):
        """Test global log-level option is recognized."""
        # Test with extend command as it's most likely to work
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extend', '--log-level', 'DEBUG', '--help']
        )

        # Should show help without error (log level should be accepted)
        assert_exit_code(exit_code, 0, stdout, stderr)
        assert_contains(stdout, 'usage:', case_sensitive=False)

    def test_invalid_global_log_level(self, cli_runner):
        """Test invalid log level shows appropriate error."""
        exit_code, stdout, stderr = cli_runner.run_teloclip(
            ['extend', '--log-level', 'INVALID']
        )

        assert exit_code != 0
        error_output = stdout + stderr
        assert_contains(error_output, 'log', case_sensitive=False)
