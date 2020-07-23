from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='SimpleNexusCaller',
    version='1.0.0',
    author='Brad Balderson, Mikael Boden',
    author_email='brad.balderson@uqconnect.edu.au',
    packages=find_packages(),
    license='GPL-3.0',
    long_description_content_type="text/markdown",
    long_description=long_description,
    scripts=['bin/simplenexuscaller'],
    install_requires = ['numpy','pandas'],
    entry_points={
        'console_scripts': [
        'simplenexuscaller=simplenexuscaller.__main__:main',
        'SimpleNexusCaller=simplenexuscaller.__main__:main'
        ]
    },
    python_requires='>=3',
    description=("SimpleNexusCaller calls ChIP-nexus peaks based on the commonly " 
				 "provided bedGraph input format. This is performed in 3 simple "
				 "steps: 1) identification of 'signal' regions on the + and - "
				 "strands, 2) identification of TF boundaries on the + and - "
				 "strand indicated by the summit of a signal range, and 3) by "
				 "matching the TF boundaries on the + strand to the closest TF "
				 "boundary downstream on the - strand."),
    keywords='SimpleNexusCaller',
    url='https://github.com/BradBalderson/SimpleNexusCaller',
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    ],
)
