from setuptools import setup

pypi_classifiers = [
    'Programming Language :: Python :: 2.7',
    'Programming Language :: Python :: 3',
    "Development Status :: 4 - Beta",
    "Environment :: Console",
    "Operating System :: OS Independent",
    'Intended Audience :: Science/Research',
    'Natural Language :: English',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    "Topic :: Software Development :: Libraries :: Python Modules",
    'License :: OSI Approved :: MIT License',
]

desc = """Filter SAM file for soft-clipped alignments containing unassembled telomeric repeats."""

setup(name='teloclip',
      version='0.0.2',
      description=desc,
      url='https://github.com/Adamtaranto/teloclip',
      author='Adam Taranto',
      author_email='aptaranto@ucdavis.edu',
      license='MIT',
      packages=['teloclip'],
      classifiers=pypi_classifiers,
      keywords=["samfile","telomere","alignment"],
      include_package_data=True,
      zip_safe=False,
      entry_points={
        'console_scripts': [
            'teloclip=teloclip.run_self:main',
        ],
    },
    )