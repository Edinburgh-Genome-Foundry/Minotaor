<p align="center">
<img alt="Minotaor logo" title="Minotaor" src="https://raw.githubusercontent.com/Edinburgh-Genome-Foundry/Minotaor/main/images/minotaor.png" width="120">
</p>


# Minotaor

[![Build Status](https://travis-ci.org/Edinburgh-Genome-Foundry/Minotaor.svg?branch=main)](https://travis-ci.org/Edinburgh-Genome-Foundry/Minotaor)
[![Coverage Status](https://coveralls.io/repos/github/Edinburgh-Genome-Foundry/Minotaor/badge.svg?branch=main)](https://coveralls.io/github/Edinburgh-Genome-Foundry/Minotaor?branch=main)

**Work in progress**

Minotaor is an a**mino** acid sequence anno**ta**t**or** for quickly identifying common protein tags and linkers in an ORF. Additionally, it can flag peptide sequences or patterns that are known to cause problems during translation. It uses Biopython.


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

# Define a custom translator:
from dna_features_viewer import BiopythonTranslator
class MinotaorTranslator(BiopythonTranslator):
    """Custom translator.

    Color warnings in red, CDS in default color, all other features in blue.
    """

    def compute_feature_color(self, feature):
        if feature.type == "warning":
            return "red"
        elif feature.type == "CDS":
            return "#7245dc"
        else:
            return "blue"

    def compute_feature_label(self, feature):
        return feature.id

graphic_record = MinotaorTranslator().translate_record(protein_record)
ax, _ = graphic_record.plot(figure_width=10, strand_in_label_threshold=7)
graphic_record.plot_sequence(ax)
```
![Example](images/example.png)


## Versioning

Minotaor uses the [semantic versioning](https://semver.org) scheme.


## License = MIT

Minotaor is [free software](https://www.gnu.org/philosophy/free-sw.en.html), which means the users have the freedom to run, copy, distribute, study, change and improve the software.

Minotaor was written at the [Edinburgh Genome Foundry](https://edinburgh-genome-foundry.github.io/) by [Peter Vegh](https://github.com/veghp) and is released under the MIT license.
