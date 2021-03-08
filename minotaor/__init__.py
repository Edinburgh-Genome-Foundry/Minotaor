from .minotaor import (
    seq_dataset,
    MinotaorTranslator,
    annotate_record,
    create_and_annotate_record,
    convert_dna_to_aa_pattern,
    make_regex_from_dna,
    convert_prosite_to_regex,
    tokenize_simple_regex,
    convert_regex_to_prosite,
    add_scanprosite_results,
)

from .version import __version__
