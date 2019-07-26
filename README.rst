==============
**randMAG**
==============

- Eric Hugoson (eric@hugoson.org / hugoson@evolbio.mpg.de / @EricHugo)


Introduction
--------------
With an ever increasing rate of sequencing comes challenges in analysing the data.
In particular, single-amplified genomes and metagenome assembled genomes (SAGs & MAGs) may
prove challenging due to their propensity to be incomplete or contaminated
with contigs from other genomes. Many tools exist to tackle the challenges at
various stages of the analysis, but none are perfect. In order to devise new
tools methods for dealing with SAGs and MAGs the ability to simulate such data
sets is paramount.

randMAG is a slim software aimed at producing artificial MAGs or SAGs from existent
complete genomes. randMAG produces randomised MAGs based on whole genomes provided, a
distribution of contig lengths given, as well as desired completeness- and
contamination levels.

Packaged with randMAG is a file containing the lengths of 192243 contigs of 2284
bacterial MAGs [1] from the Tara Oceans metagenomic survey.

#. [TullyEtAl_2018]_

Dependencies
--------------
Python (>=3.6)

Python libraries
^^^^^^^^^^^^^^^^^^^
If built from the PyPI package these will be installed automatically, otherwise can
easily be installed using ``pip`` (`Install pip <https://pip.pypa.io/en/stable/installing/>`_).

Required
""""""""""""""""""
- Biopython (>= 1.70)
- Numpy (>= 1.13.1)
- Matplotlib (>= 2.0.2)
- SciPy (>= 1.1.0)
- seaborn (>= 0.9.0)
- pandas (>= 0.23.4)

Positional arguments
^^^^^^^^^^^^^^^^^^^^^^^
::

    genome_tab            Genomes in .fna in list format
    distribution          Set of lengths to base the distribution on

``genome_tab`` should be a simple file with a list of genomes (nucleotide fasta files)
that will be used to built the artificial MAGs, e.g.: ::

    GCF_000007265.1_ASM726v1_genomic.fna
    GCF_000007325.1_ASM732v1_genomic.fna
    ...

``distribution`` should also be a simple file containing a list of integers representing
a distribution of contig sizes by which the artificial contigs will be split, e.g.::

    23967
    15981
    149609
    ...

Optionally, ``distribution`` can be filled with ``Tara_bact`` to use the set of contig
lengths from Tara Oceans described in the introduction.


Optional arguments
^^^^^^^^^^^^^^^^^^^^^^^^
  -h, --help            show this help message and exit
  -c COMPLETENESS, --completeness COMPLETENESS
                        Range of completeness levels to be produced.
                        Default=1.0
  -r CONTAMINATION, --contamination CONTAMINATION
                        Range of contamination to be included in produced
                        MAGs. Default=1.0
  -n NUM, --num NUM     Number of randomised MAGs to produce. Produces one
                        randomised per provided genome by default.


Examples
^^^^^^^^^^^^^^^^^^^^^^^^
Initially, create a file with a list of the genomic fasta files which are to be the basis for the artificial MAGs::

    $ head -n5 reference_fna.list 
    /nobackup/genomic_sources/genbank_reference_bacteria_fna_20181128/ncbi-genomes-2018-11-28/GCF_000005845.2_ASM584v2_genomic.fna
    /nobackup/genomic_sources/genbank_reference_bacteria_fna_20181128/ncbi-genomes-2018-11-28/GCF_000006745.1_ASM674v1_genomic.fna
    /nobackup/genomic_sources/genbank_reference_bacteria_fna_20181128/ncbi-genomes-2018-11-28/GCF_000006765.1_ASM676v1_genomic.fna
    /nobackup/genomic_sources/genbank_reference_bacteria_fna_20181128/ncbi-genomes-2018-11-28/GCF_000006785.2_ASM678v2_genomic.fna
    /nobackup/genomic_sources/genbank_reference_bacteria_fna_20181128/ncbi-genomes-2018-11-28/GCF_000006845.1_ASM684v1_genomic.fna


Then a second file which is a list of lengths representing the distribution of contig sizes::

    $ head -n5 Tara_bact
    23967
    15981
    149609
    37636
    121320

Example 1 - Split into contigs
""""""""""""""""""""""""""""""

This example splits up all .fna files in ``reference_fna.list`` into contigs according
to the distribution of contig lengths in ``Tara_bact``::

   $ randMAG reference_fna.list Tara_bact

The result is all a single .fna file created in the current working directory for
each input .fna. Also produced is ``simulated_MAGs.tab``, a tabular file containing input
file, the unique 8-char suffix assigned, completeness-, and redundancy fractions::

    $ head -n5 simulated_MAGs.tab
    GCF_000005845.2_ASM584v2_genomic	mpmwoutk	1.0	1.0
    GCF_000006745.1_ASM674v1_genomic	mmqehxfx	1.0	1.0
    GCF_000006765.1_ASM676v1_genomic	engiyxoz	1.0	1.0
    GCF_000006785.2_ASM678v2_genomic	rtvayhqu	1.0	1.0
    GCF_000006845.1_ASM684v1_genomic	bzpsqnfq	1.0	1.0

The content of the fasta files in this example remained unchanged apart from
being split into contigs.

Example 2 - Alter completeness/contamination
""""""""""""""""""""""""""""""""""""""""""""
To change the completeness and contamination of the fasta files the ``-c`` and ``-r``
arguments need to be used. The ``-n`` argument can be used to get precisely the
desired number of unqiue MAGs::

    $ randMAG reference_fna.list Tara_bact -c 0.7 -r 1.2 -n 10000

This will produce files that are at most 70% complete and at least 20% contaminated::

    $ head -n5 simulated_MAGs.tab 
    GCF_000005845.2_ASM584v2_genomic	kfcckaxy	0.6956260400391929	1.2092498368299183
    GCF_000006745.1_ASM674v1_genomic	xqnerzfy	0.6927500102156292	1.202911467089386
    GCF_000006765.1_ASM676v1_genomic	kiubfyau	0.6988059837775469	1.200510423558475
    GCF_000006785.2_ASM678v2_genomic	xxltcsmv	0.6849106013550827	1.2144702270571384
    GCF_000006845.1_ASM684v1_genomic	rcxezdoq	0.6927952822804169	1.2028640216997193

As well as 10 000 unique MAGs as requested with ``-n``::

    $ wc -l simulated_MAGs.tab
    10000 simulated_MAGs.tab

References
----------------

.. [TullyEtAl_2018] Tully,B.J. et al. (2018) The reconstruction of 2,631 draft metagenome-assembled genomes from the global oceans. Sci. Data, 5, 170203