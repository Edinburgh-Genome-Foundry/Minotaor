import pytest
from unittest import mock

import os

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import minotaor


testdata = [
    "F-G-Q.",
    "[GSTALIVN]-{PCHR}-{KND}-H-E-[LIVMFYW]-{DEHRKP}-H-{EKPC}-[LIVMFYWGSPQ].",
    "P-x(2)-R-G-[STAIV](2)-x-N-[APK]-x-[DE].",
    "[IV]-x-[IV]-[SA]-T-[NQ]-M-A-G-R-G-x-D-I-x-L.",
    "[GLES]-x-[LIVM]-x(2)-L-[KR]-[KRHNS]-x-K-x(5)-[LIVM]-x(2)-[GNKADS]-x-[DEN]-[CRG]-[GI].",
    "[DNSTAGC]-[GSTAPIMVQH]-x(2)-G-[DE]-S-G-[GS]-[SAPHV]-[LIVMFYWH]-[LIVMFYSTANQH].",
    "C-x(1,2)-C-x(5)-G-x(2)-C-x(2)-C-x(3,4)-[FYW]-x(3,15)-C.",
    "A-G-Y-G-S-T-x-T.",
    "[IVAG]-x-[KR]-x(2)-[DE]-[GDE](2)-x(1,2)-[EQHF]-x-[LIV]-x(4)-P-x-[LIVM](2)-[TACS].",
    "[SAPG]-[LIVMST]-[CS]-[STACG]-P-[STA]-R-x(2)-[LIVMFW](2)-[TAR]-G.",
    "C-{C}(6)-C-{C}(5)-C-C-x(1,3)-C-C-x(2,4)-C-x(3,10)-C.",
    "[LIVM]-x-[LIVMFYW]-E-G-x-[LSI]-L-K-[PA]-[SN].",
    "[FYWS]-[RK]-x-G-F-F-x-R.",
    "S-x-[LIVMF]-K-R-x(4)-K-D-x-[GSA]-x(2)-[LIF]-[PGS]-x-H-G-G-[LIVMF]-x-D-R-[LIVMFT]-D.",
    "[SA]-[FY]-[LIV]-L-[STN]-E-S-S-[LIVMF]-F-[LIV].",
    "[GR]-C-[IV]-G-R-[ILS]-x-W.",
    "<F-G(3)-x-x(0,1)-[GSTV]-[ST](2,3)-{PR}-{L}(2,5)-x-A>.",
]

mock_elm_tsv_path = os.path.join("tests", "data", "mock_elm.tsv")


def test_MinotaorTranslator():
    try:
        import dna_features_viewer
    except ImportError:
        with pytest.raises(Exception):
            minotaor.MinotaorTranslator()

    else:
        minotaor.MinotaorTranslator()


def test_annotate_record():
    protein = Seq("HHHHHHDLG*EDINBURGHGENQMEFQUNDRY")
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )

    protein_record = minotaor.annotate_record(protein_record)

    assert len(protein_record.features) == 5
    assert protein_record.features[0].id == "no start codon"
    assert protein_record.features[1].id == "ATG"
    assert protein_record.features[2].id == "not a stop codon"
    assert protein_record.features[3].id == "STOP"
    assert protein_record.features[4].id == "His tag"


def test_create_and_annotate_record():
    protein_record = minotaor.create_and_annotate_record("HHHHHH")

    assert type(protein_record) == SeqRecord
    assert protein_record.features[0].id == "no start codon"


def test_generate_prefix_codons():
    with pytest.raises(ValueError):
        minotaor.minotaor.generate_prefix_codons("")  # length not 1 nor 2


def test_generate_postfix_codons():
    with pytest.raises(ValueError):
        minotaor.minotaor.generate_postfix_codons("")  # length not 1 nor 2


def test_convert_dna_to_aa_pattern():
    with pytest.raises(ValueError):
        minotaor.convert_dna_to_aa_pattern("A")  # length less than 3

    expected = [
        "L",
        "[YSNIDVTPRCFGHAL][FL]",
        "[TSAP][YS*WCFL]",
        "K",
        "[SIVQETKP*RGAL][RS]",
        "[*QKE][VDEGA]",
    ]
    patterns = minotaor.convert_dna_to_aa_pattern("CTT")
    for index, pattern in enumerate(patterns):
        assert set(pattern) == set(expected[index])


def test_convert_prosite_to_regex():
    assert (
        minotaor.convert_prosite_to_regex("[AC]-x-V-x(4)-{ED}>.")
        == "[AC][^\*]V[^\*]{4}[^ED]$"
    )
    assert (
        minotaor.convert_prosite_to_regex("<F-[AC]-x-V-x(4)-{ED}.")
        == "^F[AC][^\*]V[^\*]{4}[^ED]"
    )
    assert (
        minotaor.convert_prosite_to_regex("<F-[AC]-x-V-x(4)-{ED}.")
        == "^F[AC][^\*]V[^\*]{4}[^ED]"
    )
    assert (
        minotaor.convert_prosite_to_regex("C-x-C-x(2)-[GP]-[FYW]-x(4,8)-C.")
        == "C[^\*]C[^\*]{2}[GP][FYW][^\*]{4,8}C"
    )
    assert (
        minotaor.convert_prosite_to_regex("<A-x-[ST](2)-x(0,1)-V.")
        == "^A[^\*][ST]{2}[^\*]{0,1}V"
    )
    with pytest.raises(ValueError):
        minotaor.convert_prosite_to_regex("V")  # no period
    with pytest.raises(Exception):
        minotaor.convert_prosite_to_regex("V-[D>].")  # '>' inside square brackets


def test_tokenize_simple_regex():
    regex = "^A[GV][^PR][FYW]{2}[^P]{4}[^\\*][^\\*]{8}$"
    assert minotaor.tokenize_simple_regex(regex) == [
        "^",
        "A",
        "[GV]",
        "[^PR]",
        "[FYW]",
        "{2}",
        "[^P]",
        "{4}",
        "[^\\*]",
        "[^\\*]",
        "{8}",
        "$",
    ]
    with pytest.raises(Exception):
        minotaor.tokenize_simple_regex("[A}")  # incorrect regex


def test_convert_tokens_to_prosite():
    tokens = [
        "^",
        "A",
        "[GV]",
        "[^PR]",
        "[FYW]",
        "{2}",
        "[^P]",
        "{4}",
        "[^\\*]",
        "[^\\*]",
        "{8}",
        "$",
    ]
    assert (
        minotaor.minotaor.convert_tokens_to_prosite(tokens)
        == "<A-[GV]-{PR}-[FYW](2)-{P}(4)-x-x(8)>."
    )


def test_convert_regex_to_prosite():
    for prosite in testdata:
        assert prosite == minotaor.convert_regex_to_prosite(
            minotaor.convert_prosite_to_regex(prosite)
        )


def test_add_scanprosite_results():
    mock_scanprosite_result = [
        {
            "sequence_ac": "USERSEQ1",
            "start": 24,
            "stop": 37,
            "signature_ac": "PS00905",
            "level_tag": "(0)",
        }
    ]
    protein = Seq("SYYHHHHHHDYDIPTTENLYFAGDLPGLMDGAAAGGGAA")
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )
    protein_record = minotaor.add_scanprosite_results(
        protein_record, mock_scanprosite_result
    )
    assert len(protein_record.features) == 1
    assert protein_record.features[0].id == "PROSITE:PS00905"
    assert int(protein_record.features[0].location.start) == 23
    assert int(protein_record.features[0].location.end) == 37


def test_get_content():
    pass


def test_evaluate_content():
    pass


def test_add_aa_content():
    my_seq = "DDDAAADDDDDAAADAAA"  # non-overlapping and overlapping
    protein_record = SeqRecord(
        Seq(my_seq), id="example", annotations={"molecule_type": "protein"}
    )
    protein_record = minotaor.add_aa_content(
        protein_record, aa=["A", "H"], window_size=5, cutoff=0.6
    )
    assert len(protein_record.features) == 2
    assert protein_record.features[0].location.start == 1
    assert protein_record.features[0].location.end == 8
    assert protein_record.features[1].location.start == 9
    assert protein_record.features[1].location.end == 18


def test_add_interpro():
    mock_interpro = mock.Mock()
    hit = mock.Mock()
    fragment = mock.Mock()

    mock_interpro.hits = [hit]
    hit.attributes = {"Hit type": "phobius"}
    hit.fragments = [fragment]
    fragment.query_start = 2
    fragment.query_end = 5
    fragment.hit_id = "Mock_ID"
    fragment.hit_description = "Mock_description"

    protein = Seq("HHHHHHDLGEDINBURGHGENQMEFQUNDRY")
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )

    protein_record = minotaor.add_interpro(protein_record, mock_interpro)
    assert len(protein_record.features) == 1
    assert protein_record.features[0].id == "phobius: Mock_ID Mock_description"


def test_add_elm_tsv():
    # Also tests add_elm()
    protein = Seq("HHHHHHDLGEDINBURGHGENQMEFQUNDRY")
    seqrecord = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )
    seqrecord = minotaor.add_elm_tsv(seqrecord, elm_tsv=mock_elm_tsv_path)
    assert len(seqrecord.features) == 2
    assert seqrecord.features[0].id == "ELM:DEG_APCC_DBOX_1"
    assert seqrecord.features[0].location.start == 6
    assert seqrecord.features[0].location.end == 14
