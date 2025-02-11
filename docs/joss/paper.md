---
title: 'Sugar: A Python framework for bioinformatics'
tags:
  - Python
  - bioinformatics
  - genomics
  - annotations
authors:
  - name: Tom Eulenfeld
    orcid: 0000-0002-8378-559X
    affiliation: 1
affiliations:
 - name: Bioinformatics/High-Throughput Analysis, Friedrich Schiller University Jena, Germany
   index: 1
date: 30 January 2025
bibliography: paper.bib

---

# Summary

Modern bioinformatics requires the use of a variety of tools [e.g. BLAST for homology detection, @blast] and databases [e.g. GenBank, @genbank]. Custom scripts often need to use the results of these tools, which are available in different formats.

``Sugar`` is a Python framework for bioinformatics and aims to facilitate rapid application development.
The package allows to read and write various sequence and annotation formats, i.e. FASTA, GenBank, Stockholm, GFF, GTF, BLAST and others. Since ``sugar`` uses a plugin system for reading and writing, new file formats can be added not only within the ``sugar`` package, but also within other packages, allowing for low barrier inclusion of new formats.
``Sugar`` includes classes for representing DNA/RNA sequences and annotations. The main functionality is exposed through methods of these classes and is therefore readily available.

``Sugar`` can be used by researchers and students of bioinformatics alike. The package has already been used as a library in the AnchoRNA package [@anchorna] and in several scripts that form the basis of the @eve_in_bat publication.
``Sugar`` can be installed from PyPI.
Online documentation and tutorials are available on the GitHub project site.

# Statement of need

Other well-known frameworks for bioinformatics are Biopython [@biopython] and Biotite [@biotite].
Launched in 2000, Biopython contains a large collection of freely available tools. In contrast, ``sugar`` tries to focus on the basics of IO as well as sequence and annotation (resp. feature) manipulation. In Table \ref{code}, we compare the code for reading a FASTA sequence file using ``sugar`` and Biopython, respectively. In addition, the ``sugar`` code example demonstrates a typical bioinformatics file manipulation task, which is not easily possible with Biopython alone due to its lack of a GFF writer. The excellent Biotite package also handles DNA/RNA sequences and has an additional focus on protein structures.
``Sugar`` should not be seen as an alternative to these other packages, but rather as a complement. Therefore, adapters for these packages are provided, which can easily convert ``sugar`` sequence objects into the corresponding objects of the Biopython and Biotite packages, and vice versa.
In addition, ``sugar`` annotation objects can be converted
to corresponding objects of the DNA features viewer package [@ftsviewer] for easy plotting.

------                                      ------
`from sugar import read, read_fts`{.python} `from Bio import SeqIO`{.python}
`seqs = read('seqs.fasta')`{.python}        `seqs = list(SeqIO.parse('seqs.fasta', 'fasta'))`{.python}
`fts = read_fts('hits.blastn')`{.python}
`seqs.fts = fts`{.python}
`seqs.write('seqs_annotated.gff')`{.python}
-----                                       ------
Table: Comparison of code to read a FASTA file using ``sugar`` (left) and Biopython (right).
The code example using ``sugar`` also demonstrates reading of a BLAST result file, attaching the hit features to the sequences while discarding features belonging to different sequences, and writing the sequences with the corresponding features to a GFF file including a FASTA directive. \label{code}

# Acknowledgements

The use of a plugin system for reading and writing was inspired by ObsPy [@obspy], a Python framework for seismology.
This project was funded by  by the Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under Germany's Excellence Strategy -- EXC 2051 -- Project ID 390713860.


# References
