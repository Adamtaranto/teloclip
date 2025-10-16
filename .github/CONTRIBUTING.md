# Contributing to teloclip

Thank you for your interest in contributing to teloclip! This document provides guidelines and instructions for contributing to this project.

## Table of Contents
- [Code of Conduct](#code-of-conduct)
- [How to Contribute](#how-to-contribute)
  - [Reporting Bugs](#reporting-bugs)
  - [Suggesting Features](#suggesting-features)
  - [Contributing Code](#contributing-code)
- [Development Setup](#development-setup)
- [Style Guidelines](#style-guidelines)
- [Testing](#testing)
- [Pull Request Process](#pull-request-process)
- [Questions and Contact](#questions-and-contact)

## Code of Conduct

Our project is committed to providing a welcoming and inclusive experience for everyone. We expect all participants to adhere to our Code of Conduct in all project spaces.

## How to Contribute

### Reporting Bugs

We use GitHub issues to track bugs. Report a bug by [opening a new issue](https://github.com/Adamtaranto/teloclip/issues/new). Be sure to include:

- A clear, descriptive title
- Steps to reproduce the issue
- Expected behavior and actual behavior
- Any relevant logs, error messages, or screenshots
- The version of teloclip you're using
- Your operating system and Python version

### Suggesting Features

We welcome feature suggestions! To suggest a feature:

1. Check the [issues page](https://github.com/Adamtaranto/teloclip/issues) to see if the feature has already been suggested
2. If not, [open a new issue](https://github.com/Adamtaranto/teloclip/issues/new) with the label "enhancement"
3. Clearly describe the feature and its potential benefits
4. If possible, outline how it might be implemented

### Contributing Code

1. Fork the repository
2. Create a new branch for your feature or bugfix
3. Make your changes
4. Run tests to ensure they pass
5. Submit a pull request

## Development Setup

To set up your development environment:

1. Clone your fork of the repository

   ```bash
   git clone https://github.com/adamtaranto/teloclip.git
   cd teloclip
  ```

3. Create and activate a virtual environment

  ```bash
  conda env create -f environment.yml
  conda activate teloclip
  ```

3. Install development dependencies

  ```bash
  pip install -e ".[dev]"
  ```

## Style Guidelines

We follow PEP 8 style guidelines and use NumPy-style docstrings.

Key points:

- Maximum line length of 88 characters
- Use type hints
- Document functions and classes with NumPy-style docstrings
- Use descriptive variable names
- Format code with Black
- Sort imports with isort

## Testing

Before submitting a pull request:

1. Ensure all tests pass

```bash
pytest tests/
```

2. For new functionality, add appropriate tests
3. Ensure test coverage remains high

## Pull Request Process

1. Update documentation if needed
2. Update the README.md if needed
3. Update the CHANGELOG.md with a description of your changes
4. Ensure all tests pass and code quality checks succeed
5. Submit your pull request with a clear description of the changes
6. Address any feedback from maintainers

## Questions and Contact

If you have questions about contributing, please [open an issue](https://github.com/Adamtaranto/teloclip/issues) or contact the maintainer directly.

Thank you for contributing to teloclip!
