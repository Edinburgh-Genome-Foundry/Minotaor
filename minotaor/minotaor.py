import os
import re
import pandas

from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

DATA_DIR = os.path.join(os.path.dirname(os.path.realpath(__file__)), "data")
SEQ_DATA = pandas.read_csv(os.path.join(DATA_DIR, "seq.csv"))

try:
    from dna_features_viewer import BiopythonTranslator
except ImportError:

    class MinotaorTranslator:
        """Please install dna_features_viewer to use this class."""

        def __init__(self):
            raise Exception("Please install dna_features_viewer to use this class.")


else:

    class MinotaorTranslator(BiopythonTranslator):
        """Custom translator.

        Color warnings in red, CDS in default color, all other features in blue.
        """

        def compute_feature_color(self, feature):
            mino_class = "mino_class"
            if feature.qualifiers[mino_class] == "error":
                return "red"
            elif feature.qualifiers[mino_class] == "warning":
                return "yellow"
            elif feature.qualifiers[mino_class] == "tag":
                return "tab:blue"
            elif feature.qualifiers[mino_class] == "linker":
                return "tab:cyan"
            else:
                return "#7245dc"  # default dna_features_viewer color

        def compute_feature_label(self, feature):
            try:
                label = feature.qualifiers["label"]
            except Exception:
                label = feature.id

            return label


def annotate_record(seqrecord, seq_dataset=None):
    """Annotate a record with entries of a reference sequence dataset.


    Note that the search is case sensitive.

    **Parameters**

    **seqrecord**
    > SeqRecord to annotate.

    **seq_dataset**
    > A minotaor sequence dataset (`pandas.DataFrame`). Default uses the built-in data.
    The `sequence` and `name` columns are used for search and naming of the motifs.
    If there is a `class` column, then it is used for the SeqFeature's mino_class
    qualifier, which is used to determine the color during plotting. If the dataframe
    has a `description` column, then its entry is added as a `note` qualifier.
    """
    if seq_dataset is None:
        seq_dataset = SEQ_DATA
    # FLAG NO START: M
    if str(seqrecord.seq)[0] != "M":
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(0, 1),
                type="misc_feature",
                id="no start codon",
                qualifiers={"label": "no start codon", "mino_class": "warning"},
            )
        )
    # FLAG INTERNAL START: M (ATG)
    atg_positions = [i for i, letter in enumerate(str(seqrecord.seq)) if letter == "M"]
    for position in atg_positions:
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(position, position + 1),
                type="misc_feature",
                id="ATG",
                qualifiers={"label": "ATG", "mino_class": "warning"},
            )
        )
    # FLAG NO END: *
    if str(seqrecord.seq)[-1] != "*":
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(len(seqrecord) - 1, len(seqrecord)),
                type="misc_feature",
                id="not a stop codon",
                qualifiers={"label": "not a stop codon", "mino_class": "warning"},
            )
        )
    # FLAG STOP CODONS: *
    stop_positions = [i for i, letter in enumerate(str(seqrecord.seq)) if letter == "*"]
    for position in stop_positions:
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(position, position + 1),
                type="misc_feature",
                id="STOP",
                qualifiers={"label": "STOP", "mino_class": "error"},
            )
        )
    # ANNOTATE PATTERNS
    if "class" in seq_dataset.columns:
        has_mino_class = True
        class_list = seq_dataset["class"].to_list()
    else:
        has_mino_class = False
        mino_class = "default"
    if "description" in seq_dataset.columns:
        has_description = True
        description_list = seq_dataset["description"].to_list()
    else:
        has_description = False
    patterns = seq_dataset["sequence"].to_list()
    names = seq_dataset["name"].to_list()
    for index, pattern in enumerate(patterns):
        name = names[index]
        if has_mino_class:
            mino_class = class_list[index]
        matches = {m.start(): m.end() for m in re.finditer(pattern, str(seqrecord.seq))}

        if has_description:
            description = description_list[index]
            note = name + ". Description: " + str(description)
        else:
            note = name
        qualifier = {"label": name, "mino_class": mino_class, "note": note}

        for start, end in matches.items():
            seqrecord.features.append(
                SeqFeature(
                    FeatureLocation(start, end),
                    type="misc_feature",
                    id=name,
                    qualifiers=qualifier,
                )
            )

    return seqrecord


def create_and_annotate_record(sequence, seq_dataset=None):
    """Create a SeqRecord from an amino acid sequence string.


    **Parameters**

    **sequence**
    > Sequence (`str`).
    """
    if seq_dataset is None:
        seq_dataset = SEQ_DATA
    protein = Seq(sequence)
    protein_record = SeqRecord(
        protein, id="example", annotations={"molecule_type": "protein"}
    )
    protein_record = annotate_record(protein_record)

    return protein_record


def convert_dna_to_aa_pattern(dna):
    """Convert a DNA string to a list of patterns representing its translations.


    **Parameters**

    **dna**
    > DNA (`str` of `ATCG`)"""
    if len(dna) < 3:
        raise ValueError("Minimum DNA length is 3")

    patterns = []

    for frame in [0, 1, 2]:
        aa_dna = dna[frame:]
        prefix = dna[:frame]
        modulo = len(aa_dna) % 3  # codon length is 3
        if modulo != 0:
            postfix = aa_dna[-modulo:]
            aa_dna = aa_dna[:-modulo]
        else:
            postfix = ""

        regex = make_regex_from_dna(aa_dna, prefix, postfix)
        patterns += [regex]

    dna_reverse_complement = str(Seq(dna).reverse_complement())
    for frame in [0, 1, 2]:
        aa_dna = dna_reverse_complement[frame:]
        prefix = dna_reverse_complement[:frame]
        modulo = len(aa_dna) % 3  # codon length is 3
        if modulo != 0:
            postfix = aa_dna[-modulo:]
            aa_dna = aa_dna[:-modulo]
        else:
            postfix = ""

        regex = make_regex_from_dna(aa_dna, prefix, postfix)
        patterns += [regex]

    return patterns


def make_regex_from_dna(dna, prefix, postfix):
    """Convert three DNA strings into a regex.

    The first DNA string (`dna`) must be divisible by 3, the length of the
     second (`prefix`) and third (`postfix`) must be 1 or 2.


    **Parameters**

    **dna**
    > DNA (`str` of `ATCG`).

    **prefix**
    > DNA (`str` of `ATCG`).

    **postfix**
    > DNA (`str` of `ATCG`).
    """

    aa = str(Seq(dna).translate())
    prefix_regex = create_prefix_regex(prefix)
    postfix_regex = create_postfix_regex(postfix)
    regex = prefix_regex + aa + postfix_regex

    return regex


def create_prefix_regex(prefix):
    if prefix:
        prefix_codons = generate_prefix_codons(prefix)
        translated_prefixes = []
        for codon in prefix_codons:
            aa = str(Seq(codon).translate())
            translated_prefixes += [aa]
        translated_prefixes = list(set(translated_prefixes))  # remove duplicates
        prefix_regex = "".join(translated_prefixes)
        prefix_regex = "[" + prefix_regex + "]"  # match 1
    else:
        prefix_regex = ""
    return prefix_regex


def create_postfix_regex(postfix):
    if postfix:
        postfix_codons = generate_postfix_codons(postfix)
        translated_postfixes = []
        for codon in postfix_codons:
            aa = str(Seq(codon).translate())
            translated_postfixes += [aa]
        translated_postfixes = list(set(translated_postfixes))  # remove duplicates
        postfix_regex = "".join(translated_postfixes)
        postfix_regex = "[" + postfix_regex + "]"  # match 1
    else:
        postfix_regex = ""

    return postfix_regex


def generate_prefix_codons(prefix):
    codons = []
    if len(prefix) == 1:
        for first_letter in ["A", "T", "C", "G"]:
            for second_letter in ["A", "T", "C", "G"]:
                codon = first_letter + second_letter + prefix
                codons += [codon]
    elif len(prefix) == 2:
        for first_letter in ["A", "T", "C", "G"]:
            codon = first_letter + prefix
            codons += [codon]
    else:
        raise ValueError("Length of prefix must be 1 or 2")

    return codons


def generate_postfix_codons(postfix):
    codons = []
    if len(postfix) == 1:
        for second_letter in ["A", "T", "C", "G"]:
            for third_letter in ["A", "T", "C", "G"]:
                codon = postfix + second_letter + third_letter
                codons += [codon]
    elif len(postfix) == 2:
        for third_letter in ["A", "T", "C", "G"]:
            codon = postfix + third_letter
            codons += [codon]
    else:
        raise ValueError("Length of postfix must be 1 or 2")

    return codons


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


def add_scanprosite_results(seqrecord, scanprosite_record):
    for entry in scanprosite_record:
        start_index = entry["start"] - 1  # convert to Python indexing
        stop_index = entry["stop"]  # prosite range inclusive, Python not
        name = "PROSITE:" + entry["signature_ac"]  # key for the prosite ID
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(start_index, stop_index),
                type="misc_feature",
                id=name,
                qualifiers={"label": name, "mino_class": "default"},
            )
        )
    return seqrecord


def get_content(sequence, aa, window_size):
    """Compute proportion of selected amino acids in string."""
    proportions = []
    for index in range(0, (len(sequence) - window_size + 1)):
        subseq = sequence[index : index + window_size]
        occurrences = 0
        for letter in aa:
            occurrences += subseq.count(letter)
        proportion = occurrences / window_size
        proportions += [proportion]

    return proportions


def evaluate_content(sequence, aa, window_size, cutoff):
    """Compute global or local content of selected amino acids."""
    proportions = get_content(sequence, aa, window_size)

    positions = []
    for index, proportion in enumerate(proportions):
        if proportion >= cutoff:
            positions += [index]  # add to breaches

    # sum positions:
    ranges = {}
    ranges[positions[0]] = positions[0] + window_size  # initialize
    previous = positions[0]
    for position in positions[1:]:
        if ranges[previous] >= position:
            new_end = position + window_size
            ranges[previous] = new_end
            # key remains `previous`
        else:
            end = position + window_size
            ranges[position] = end
            previous = position

    return positions, ranges


def add_aa_content(seqrecord, aa, window_size, cutoff, name=None):
    """Compute and annotate global or local content of selected amino acids.


    **Parameters**

    **seqrecord**
    > The amino acid SeqRecord to annotate.

    **aa**
    > List of amino acids to search for (`list`).

    **window_size**
    > The search window size (`int`).

    **cutoff**
    > Annotate section with at least this proportion (between 0 and 1).

    **name**
    > Annotation label (`str`). Default: `>=#% X/Y/Z`.
    """
    if name is None:
        name = ">=" + str(int(cutoff * 100)) + "% " + "/".join(aa)
    sequence = str(seqrecord.seq)
    positions, ranges = evaluate_content(
        sequence, aa=aa, window_size=window_size, cutoff=cutoff
    )
    for start, stop in ranges.items():
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(start, stop),
                type="misc_feature",
                id=name,
                qualifiers={"label": name, "mino_class": "warning"},
            )
        )

    return seqrecord


def add_interpro(seqrecord, interpro, hit_types=None, include_description=True):
    """Annotate SeqRecord with InterPro results.


    **Parameters**

    **seqrecord**
    > `SeqRecord` to annotate.

    **interpro**
    > `QueryResult` object output of `Bio.SearchIO.read(handle, "interproscan-xml")`.

    **hit_types**
    > The InterProScan hit types to filter for (`list`). Default includes all.

    **include_description**
    > If True, includes description in the label, otherwise only in the `note`
    qualifier of the SeqRecord.
    """
    for hit in interpro.hits:
        if hit_types is None or hit.attributes["Hit type"] in hit_types:
            for fragment in hit.fragments:
                start = fragment.query_start
                end = fragment.query_end
                identifier = "%s: %s" % (hit.attributes["Hit type"], fragment.hit_id)
                if include_description:
                    identifier = identifier + " " + fragment.hit_description
                qualifier = {
                    "note": fragment.hit_description,
                    "label": identifier,
                    "mino_class": "default",
                }
                seqrecord.features.append(
                    SeqFeature(
                        FeatureLocation(start, end),
                        type="interpro",
                        id=identifier,
                        qualifiers=qualifier,
                    )
                )

    return seqrecord


def add_elm_tsv(seqrecord, elm_tsv):
    """Annotate SeqRecord with Eukaryotic Linear Motif (ELM) database search results.


    **Parameters**

    **seqrecord**
    > `SeqRecord` to annotate.

    **elm_tsv**
    > Path to TSV file of results (`str`).
    """
    elm = pandas.read_csv(elm_tsv, sep="\t")

    return add_elm(seqrecord, elm)


def add_elm(seqrecord, elm):
    """Annotate SeqRecord with Eukaryotic Linear Motif (ELM) database search results.


    **Parameters**

    **seqrecord**
    > `SeqRecord` to annotate.

    **elm**
    > Dataframe of results (`pandas.DataFrame`).
    """
    for row in elm.itertuples(index=True, name="Pandas"):
        seqrecord.features.append(
            SeqFeature(
                FeatureLocation(row.start - 1, row.stop - 1),  # indexing starts from 1
                type="misc_feature",
                id="ELM:" + row.elm_identifier,
                qualifiers={
                    "label": "ELM:" + row.elm_identifier,
                    "mino_class": "default",
                },
            )
        )

    return seqrecord
