from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import minotaor


def test_annotate_record():
    protein = Seq("HHHHHHDLG*EDINBURGHGENQMEFQUNDRY")
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )

    protein_record = minotaor.annotate_record(protein_record)

    assert len(protein_record.features) == 5
    assert protein_record.features[0].id == "no start codon"
    assert protein_record.features[1].id == "not a stop codon"
    assert protein_record.features[2].id == "STOP"
    assert protein_record.features[3].id == "6xHis"


def test_create_and_annotate_record():
    protein_record = minotaor.create_and_annotate_record("HHHHHH")

    assert type(protein_record) == SeqRecord
    assert protein_record.features[0].id == "no start codon"
