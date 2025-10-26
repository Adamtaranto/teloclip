#!/usr/bin/env python3
"""
Integration test validation script for teloclip extend command.

This script validates that the complete integration test framework is working
properly by running the actual tests and verifying the synthetic data generation
pipeline produces expected results.

Usage:
    python test_validate_integration.py [--quick] [--verbose]

Options:
    --quick: Run only basic validation tests
    --verbose: Show detailed output from all test runs
"""

import argparse
from pathlib import Path
import subprocess
import sys


def run_command(cmd: list, description: str, verbose: bool = False) -> bool:
    """Run a command and return success status."""
    if verbose:
        print(f'\n=== {description} ===')
        print(f'Running: {" ".join(cmd)}')

    try:
        result = subprocess.run(cmd, capture_output=not verbose, text=True, timeout=120)

        if result.returncode == 0:
            if verbose:
                print(f'✅ {description} - SUCCESS')
                if result.stdout:
                    print(f'Output:\n{result.stdout}')
            else:
                print(f'✅ {description}')
            return True
        else:
            print(f'❌ {description} - FAILED')
            if result.stderr:
                print(f'Error: {result.stderr}')
            return False

    except subprocess.TimeoutExpired:
        print(f'❌ {description} - TIMEOUT')
        return False
    except Exception as e:
        print(f'❌ {description} - ERROR: {e}')
        return False


def check_test_data_exists() -> bool:
    """Check if test data files exist."""
    test_data_dir = Path('tests/integration/test_data')
    required_files = [
        'synthetic_contigs.fasta',
        'synthetic_contigs.fasta.fai',
        'synthetic_alignments.sam',
        'synthetic_alignments_sorted.bam',
        'synthetic_alignments_sorted.bam.bai',
    ]

    print('\n📁 Checking test data files...')
    all_exist = True

    for file_name in required_files:
        file_path = test_data_dir / file_name
        if file_path.exists():
            size = file_path.stat().st_size
            print(f'  ✅ {file_name} ({size} bytes)')
        else:
            print(f'  ❌ {file_name} - MISSING')
            all_exist = False

    return all_exist


def regenerate_test_data(verbose: bool = False) -> bool:
    """Regenerate test data from scratch."""
    print('\n🔄 Regenerating test data...')

    # Change to tests/integration directory
    test_dir = Path('tests/integration')

    # Generate contigs
    success = run_command(
        ['python', 'test_generate_data.py'], 'Generate synthetic contigs', verbose
    )
    if not success:
        return False

    # Generate alignments
    success = run_command(
        ['python', 'test_generate_alignments.py'],
        'Generate synthetic alignments',
        verbose,
    )
    if not success:
        return False

    # Create BAM files
    test_data_dir = test_dir / 'test_data'

    success = run_command(
        ['samtools', 'view', '-bS', str(test_data_dir / 'synthetic_alignments.sam')],
        'Convert SAM to BAM',
        verbose,
    )
    if not success:
        return False

    # Note: The above command writes to stdout, need to redirect
    # Let's use a different approach
    cmd = [
        'bash',
        '-c',
        f'cd {test_data_dir} && samtools view -bS synthetic_alignments.sam > synthetic_alignments.bam',
    ]
    success = run_command(cmd, 'Convert SAM to BAM (redirected)', verbose)
    if not success:
        return False

    success = run_command(
        [
            'samtools',
            'sort',
            str(test_data_dir / 'synthetic_alignments.bam'),
            '-o',
            str(test_data_dir / 'synthetic_alignments_sorted.bam'),
        ],
        'Sort BAM file',
        verbose,
    )
    if not success:
        return False

    success = run_command(
        ['samtools', 'index', str(test_data_dir / 'synthetic_alignments_sorted.bam')],
        'Index sorted BAM',
        verbose,
    )
    if not success:
        return False

    return True


def validate_teloclip_functionality(verbose: bool = False) -> bool:
    """Test basic teloclip extend functionality."""
    print('\n🧪 Validating teloclip extend functionality...')

    test_data_dir = Path('tests/integration/test_data')
    fasta_file = test_data_dir / 'synthetic_contigs.fasta'
    bam_file = test_data_dir / 'synthetic_alignments_sorted.bam'

    # Run teloclip extend with custom prefix to avoid conflicts
    success = run_command(
        [
            'teloclip',
            'extend',
            str(bam_file),
            str(fasta_file),
            '--prefix',
            'validation_test',
            '--dry-run',  # Don't write output files
        ],
        'Run teloclip extend (dry run)',
        verbose,
    )

    return success


def run_integration_tests(quick: bool = False, verbose: bool = False) -> bool:
    """Run the pytest integration tests."""
    print('\n🔬 Running integration tests...')

    if quick:
        # Run just the basic functionality test
        cmd = [
            'python',
            '-m',
            'pytest',
            'tests/integration/test_teloclip_extend_integration.py::TestBasicFunctionality::test_basic_extend_runs_successfully',
            '-v',
        ]
        return run_command(cmd, 'Basic integration test', verbose)
    else:
        # Run all integration tests
        cmd = [
            'python',
            '-m',
            'pytest',
            'tests/integration/test_teloclip_extend_integration.py',
            '-v',
        ]
        return run_command(cmd, 'Full integration test suite', verbose)


def main():
    """Main validation function."""
    parser = argparse.ArgumentParser(
        description='Validate teloclip integration test framework'
    )
    parser.add_argument(
        '--quick', action='store_true', help='Run only basic validation tests'
    )
    parser.add_argument(
        '--verbose', action='store_true', help='Show detailed output from all test runs'
    )
    parser.add_argument(
        '--regenerate', action='store_true', help='Force regeneration of test data'
    )

    args = parser.parse_args()

    print('🔍 Teloclip Integration Test Validation')
    print('=' * 50)

    # Check if we're in the right directory
    if not Path('tests/integration').exists():
        print('❌ Error: Must be run from teloclip project root directory')
        return 1

    # Change to project root if needed
    if Path.cwd().name == 'integration':
        import os

        os.chdir('../..')

    success_count = 0
    total_tests = 0

    # Check test data exists or regenerate if requested
    if args.regenerate or not check_test_data_exists():
        total_tests += 1
        if regenerate_test_data(args.verbose):
            success_count += 1
        else:
            print('\n❌ Test data generation failed - cannot proceed')
            return 1

    # Validate teloclip functionality
    total_tests += 1
    if validate_teloclip_functionality(args.verbose):
        success_count += 1

    # Run integration tests
    total_tests += 1
    if run_integration_tests(args.quick, args.verbose):
        success_count += 1

    # Final summary
    print('\n' + '=' * 50)
    print(f'📊 Validation Summary: {success_count}/{total_tests} tests passed')

    if success_count == total_tests:
        print('🎉 All validation tests passed! Integration test framework is ready.')
        print('\nTo run the full test suite manually:')
        print('  cd tests/integration')
        print('  python -m pytest test_teloclip_extend_integration.py -v')
        return 0
    else:
        print('💥 Some validation tests failed. Please check the output above.')
        return 1


if __name__ == '__main__':
    sys.exit(main())
