SimpleNexusCaller
==================

**SimpleNexusCaller calls ChIP-nexus peaks** based on the commonly provided bedGraph
input format. This is performed in 3 simple steps: 1) identification of 'signal'
regions on the + and - strands, 2) identification of TF boundaries on the + and - 
strand indicated by the summit of a signal range, and 3) by matching the 
TF boundaries on the + strand to the closest TF boundary downstream on the - 
strand. **See designNotes.txt to better understand the implimentation.**

Install
-------

- [Python3.x](https://www.python.org/getit/) with the following packages:
- Numpy
- Pandas
    
To install from source:

    git clone https://github.com/BradBalderson/SimpleNexusCaller.git
    cd SimpleNexusCaller
    python3 setup.py install

Usage
-----

In the command line, type in **'simplenexuscaller -h '** for detailed usage.

    $ simplexnexuscaller -h
    
    usage: simplenexuscaller [-h] -i INPUT INPUT [-c CUTOFF] [-f FALSEINROWUPPER]
                         [-n NINROWCUTOFF] [-d DISTLIMIT] [-m MAXWIDTH]
                         [-o OUTPUT]

    Takes ChIP-nexus data in bedGraph format for + and - strand and performs fast
    and simple peak calling.
    
    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT INPUT, --input INPUT INPUT
                            ChIP-nexus bedGraph files, where each column is [chr,
                            start, end, count], with no column headers. These
                            files must be in the order of counts on the + or -
                            strand. Inputs are separated by a single space. Each
                            BedGraph contains positions per base where counts
                            represent the number of 5' reads mapping to that
                            position. Where no reads mapped, position refers to
                            interval where no reads mapped.
      -c CUTOFF, --cutoff CUTOFF
                            Cutoff number of counts above which the
                            positionconsidered as having signal.
      -f FALSEINROWUPPER, --falseInRowUpper FALSEINROWUPPER
                            No. of no signal positions (count>cutoff) in row
                            before terminate extension of signal region.
      -n NINROWCUTOFF, --nInRowCutoff NINROWCUTOFF
                            The minimum length of a TF edge signal for thesignal
                            to be called as a true signal.
      -d DISTLIMIT, --distLimit DISTLIMIT
                            Minimum distance between signal range on the
                            samestrand for them to be considered the same or
                            differentTF binding signal edges.
      -m MAXWIDTH, --maxWidth MAXWIDTH
                            Maximum width a peak is allowed to be.
      -o OUTPUT, --output OUTPUT
                            Output filename prefix. Automatically adds .bed



Example
------
    $ simplenexuscaller -i posCounts.bedGraph negCounts.bedGraph -o output_prefix   

Output
------

Output will be in standard bed file format:

- **output_prefix.bed**: The called peaks. 

output_prefix.bed file has 3 columns. See the toy example below.

|chr |start|end  |
|----|-----|-----|
|chr1|9118 |10409|

Citation
--------

Contact
-------

Authors: Brad Balderson, Mikael Boden

Contact:  brad.balderson@uqconnect.edu.au
