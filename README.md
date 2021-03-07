<p align="center">
<img alt="Minotaor logo" title="Minotaor" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Minotaor/main/images/minotaor.png" width="120">
</p>


# Minotaor

[![Build Status](https://travis-ci.org/Edinburgh-Genome-Foundry/Minotaor.svg?branch=main)](https://travis-ci.org/Edinburgh-Genome-Foundry/Minotaor)
[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/Minotaor/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/Minotaor?branch=main)

**Work in progress**

Minotaor is an a**mino** acid sequence anno**ta**t**or** for quickly identifying common
protein tags and linkers in an ORF. Additionally, it can flag peptide sequences or
patterns that are known to cause problems during translation. It uses Biopython.


## Background

In the PROSITE nomenclature, a sequence *motif* is a description of the occurrence of 
amino acids (signature, fingerprint), and can be either a pattern or a profile.
A *pattern* is a qualitative description of a motif in a regular expression-like syntax.
A *profile* (or weight matrix) is a table of position-specific amino acid weights and
gap costs.


## Install

```
pip install minotaor
```


## Usage

```python
import minotaor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

protein = Seq("SYYHHHHHHDYDIPTTENLYFQG*EDINBURGHGENQMEFQUNDRY*")
protein_record = SeqRecord(protein, id="example", annotations={"molecule_type": "protein"})

protein_record = minotaor.annotate_record(protein_record)
```

Plotting requires DNA Features Viewer installed:
```python
graphic_record = minotaor.MinotaorTranslator().translate_record(protein_record)
ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
graphic_record.plot_sequence(ax)
```
![Example](images/example.png)

Convert between PROSITE and regex formats:
```python
regex = minotaor.convert_prosite_to_regex("<A-[GV]-{PR}-[FYW](2)-{P}(4)-x-x(8)>.")
regex
# '^A[GV][^PR][FYW]{2}[^P]{4}[^\\*][^\\*]{8}$'
minotaor.convert_regex_to_prosite(regex)
# '<A-[GV]-{PR}-[FYW](2)-{P}(4)-x-x(8)>.'
```

Given a DNA sequence, return amino acid sequences that may contain it:
```python
bsmbi_site = "CGTCTC"
print(minotaor.convert_dna_to_aa_pattern(bsmbi_site))
# ['RL', '[DATRSICYNLFPHVG]V[S]', '[SPAT]S[RPLQH]', 'ET', '[MAT*RSPKLEVQGW]R[R]', '[G*R]D[DAVEG]']
```
Returns a regex for each of the 6 translation frames.


## Versioning

Minotaor uses the [semantic versioning](https://semver.org) scheme.


## License = MIT

Minotaor is [free software](https://www.gnu.org/philosophy/free-sw.en.html), which means
the users have the freedom to run, copy, distribute, study, change and improve the software.

Minotaor was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/)
by [Peter Vegh](https://github.com/veghp) and is released under the MIT license.
