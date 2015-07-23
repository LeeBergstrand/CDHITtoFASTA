CDHITtoFASTA
============
### Extracts CD-Hit clusters which contain reference proteins and stores them in FASTA format.

CD-HIT is a widely used program for clustering biological sequences to reduce sequence redundancy and improve the performance of other sequence analyses. CD-HIT was originally developed to cluster protein sequences to create reference databases with reduced redundancy and was then extended to support clustering nucleotide sequences and comparing two datasets.

CD-Hit can also be used to cluster sequences and outputs a cluster file (extension .clstr) contains the accessions of clustered sequence from an input FASTA file.

CDHITtoFASTA uses this cluster file to filter input FASTA files by extracting sequences from the file which CD-Hit finds to cluster with reference sequences.

Command Line Interface Overview
-------------------------------
```
$ python ~/PATH_TO_DIR/CDHITtoFASTA -h
usage: CDHITtoFASTA [-h] [-i CLUSTER] [-s FASTA] [-r LIST]

Extracts CD-Hit clusters which contain reference proteins and stores them in
FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  -i CLUSTER, --cluster_file CLUSTER
                        CD-Hit cluster file which provides clustering
                        information.
  -s FASTA, --sequence_file FASTA
                        FASTA file which provides sequences to be extracted.
  -r LIST, --reference_list LIST
                        File of sequence identifiers (one per line) who's CD-
                        HIT clusters should turned into FASTA files.
```

Dependencies
------------
- **Python** 2.7, 3.4 or newer
- Python Packages:
	- **biopython** 1.63
- **CD-Hit** 4.6 or newer

### To install dependancies on UNIX operating systems:

1. Install CD-Hit via Homebrew or Linuxbrew:

    ```
    brew install cd-hit
    ```
2. Install Biopython via pip:

    ```
    pip install biopython
    ```

Licence
-------

BioMagick is open-source and released under [MIT License](http://en.wikipedia.org/wiki/MIT_License).

	The MIT License (MIT)
	
	Copyright (c) 2014 Lee Bergstrand
	
	Permission is hereby granted, free of charge, to any person obtaining a copy
	of this software and associated documentation files (the "Software"), to deal
	in the Software without restriction, including without limitation the rights
	to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
	copies of the Software, and to permit persons to whom the Software is
	furnished to do so, subject to the following conditions:
	
	The above copyright notice and this permission notice shall be included in all
	copies or substantial portions of the Software.
	
	THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
	IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
	FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
	AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
	LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
	OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
	SOFTWARE.

Authors
-------
[Lee Bergstrand](http://github.com/LeeBergstrand)
