import os
import re
import pandas

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
seq_dataset = pandas.read_csv(os.path.join(DATA_DIR, "seq.csv"))


def annotate_record(seqrecord, seq_dataset=seq_dataset):
    """Annotate a record with a reference sequence dataset.


    **Parameters**

    **seqrecord**
    > SeqRecord to annotate.

    **seq_dataset**
    > A minotaor sequence dataset (`pandas.DataFrame`).
    """
    # FLAG NO START: M
    if str(seqrecord.seq)[0] != "M":
        seqrecord.features.append(
            SeqFeature(FeatureLocation(0, 1), type="warning", id="no start codon")
        )
    # FLAG NO END: *
    if str(seqrecord.seq)[-1] != "*":
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(len(seqrecord) - 1, len(seqrecord)),
                type="warning",
                id="not a stop codon",
            )
        )
    # FLAG STOP CODONS: *
    stop_positions = [i for i, letter in enumerate(str(seqrecord.seq)) if letter == "*"]
    for position in stop_positions:
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(position, position + 1), type="warning", id="STOP"
            )
        )
    # ANNOTATE SEQUENCES
    sequences = seq_dataset.loc[seq_dataset["type"] == "seq"]["sequence"].to_list()
    names = seq_dataset.loc[seq_dataset["type"] == "seq"]["name"].to_list()
    for index, sequence in enumerate(sequences):
        len_sequence = len(sequence)
        name = names[index]
        matches = [
            m.start() for m in re.finditer(re.escape(sequence), str(seqrecord.seq))
        ]
        for match in matches:
            seqrecord.features.append(
                SeqFeature(
                    FeatureLocation(match, (match + len_sequence)), type="CDS", id=name
                )
            )
    # ANNOTATE PATTERNS
    patterns = seq_dataset.loc[seq_dataset["type"] == "pattern"]["sequence"].to_list()
    names = seq_dataset.loc[seq_dataset["type"] == "pattern"]["name"].to_list()
    for index, pattern in enumerate(patterns):
        name = names[index]
        matches = {m.start(): m.end() for m in re.finditer(pattern, str(seqrecord.seq))}
        for start, end in matches.items():
            seqrecord.features.append(
                SeqFeature(FeatureLocation(start, end), type="CDS", id=name)
            )

    return seqrecord


def create_and_annotate_record(sequence, seq_dataset=seq_dataset):
    """Create a SeqRecord from an amino acid sequence string.


    **Parameters**

    **sequence**
    > Sequence (`str`).
    """
    if seq_dataset is None:
        seq_dataset = seq_dataset
    protein = Seq(sequence)
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )
    protein_record = annotate_record(protein_record)

    return protein_record


def convert_prosite_to_regex(prosite_string):
    """Convert a PROSITE motif string to a regex string.


    **Parameters**

    **prosite_string**
    > The PROSITE string (`str`).
    """
    # Implemented with a hack: by replacing characters, instead of using a lexer.
    # See https://prosite.expasy.org/prosuser.html#conv_pa for definition.
    # Remove period that ends the pattern:
    if prosite_string[-1] == ".":
        prosite_string = prosite_string[:-1]
    else:
        raise ValueError("Invalid format: a period ('.') must end the pattern")

    # N- and C-terminal restrictions:
    if "<" in prosite_string:
        N_terminal = True
        prosite_string = prosite_string.replace("<", "")
    else:
        N_terminal = False
    if ">" in prosite_string:
        if prosite_string[-1] != ">":
            raise Exception("'>' inside square brackets is not supported yet")
        C_terminal = True
        prosite_string = prosite_string.replace(">", "")
    else:
        C_terminal = False

    tokens = prosite_string.split("-")
    regex_tokens = []
    for token in tokens:
        # Convert 'x' to regex: any amino acid, but don't match stop codons.
        token = token.replace("x", "[^\\*]")

        # Replace braces for exceptions.
        if token[0] == "{":
            token = token.replace("{", "[^")
            token = token.replace("}", "]")

        # Replace for repetition. Must come after exception replacement.
        token = token.replace("(", "{")
        token = token.replace(")", "}")

        regex_tokens += [token]

    regex = "".join(regex_tokens)
    if N_terminal:
        regex = "^" + regex
    if C_terminal:
        regex += "$"

    return regex


def convert_regex_to_prosite(regex):
    """Convert a compatible regex string to a PROSITE motif.


    **Parameters**

    **regex**
    > The regex string (`str`).
    """
    tokens = tokenize_simple_regex(regex)
    # print(tokens)
    regex = convert_tokens_to_prosite(tokens)
    return regex


def convert_tokens_to_prosite(tokens):
    # The first ^ signifies N-terminal position
    if tokens[0] == "^":
        is_N_terminal = True
        del tokens[0]
    else:
        is_N_terminal = False
    regex_tokens = []
    if tokens[-1] == "$":
        is_C_terminal = True
        del tokens[-1]
    else:
        is_C_terminal = False

    # These are in reverse order compared to convert_prosite_to_regex():
    for token in tokens:
        # Replace for repetition. Must come before exception replacement.
        token = token.replace("{", "(")
        token = token.replace("}", ")")

        # Replace braces for exceptions.
        if token[0:2] == "[^":
            token = token.replace("[^", "{")
            token = token.replace("]", "}")

        # Convert wildcard to 'x'; but keep repetition (x,y)
        if "*" in token:
            token = token.replace("{", "")
            token = token.replace("}", "")
            token = token.replace("\\", "")
            token = token.replace("*", "x")

        regex_tokens = regex_tokens + [token]

    regex = tokens[0]
    for token in regex_tokens[1:]:
        if token[0] == "(":
            regex += token
        else:
            regex = regex + "-" + token

    if is_N_terminal:
        regex = "<" + regex  # < is a prosite symbol
    if is_C_terminal:
        regex = regex + ">"  # > is a prosite symbol
    # Add period that must end the pattern:
    regex = regex + "."

    return regex


def tokenize_simple_regex(regex):
    """Lex regex into list of tokens."""
    # As the format of compatible regexes are simple, a tokenizer is implemented here,
    # instead of using an external lexer with a grammar definition.
    tokens = []
    index = 0
    token_boundaries = []  # this collects the start index of each token
    closing_brackets = {"[": "]", "{": "}", "(": ")"}
    is_group = False
    for index, character in enumerate(regex):
        if character in "[{(":
            is_group = True
            token_boundaries += [index]
            closing_bracket = closing_brackets[character]
        elif not is_group:
            token_boundaries += [index]
        if character in "]})":
            if character != closing_bracket:
                raise Exception("Regex incorrect or cannot be converted to PROSITE.")
            is_group = False

    for index, boundary in enumerate(token_boundaries):
        try:
            token = regex[boundary : token_boundaries[index + 1]]
        except IndexError:
            token = regex[boundary:]  # last token

        tokens += [token]

    return tokens
