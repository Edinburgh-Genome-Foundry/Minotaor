import pandas

from Bio import SeqIO

# Include these:
keywords = [
    "Basic",
    "Reporter",
    "Translational_Unit",
    "Composite",
    "Protein_Domain",
    "Signalling",
    "Coding",
    "Tag",
    "Device",
]

aa_sequences = []
names = []
descriptions = []
# get fasta and info on header lines from http://parts.igem.org/Registry_API
for seq_record in SeqIO.parse("All_Parts.fasta", "fasta"):
    # take only short DNA that can be translated:
    if 8 < len(seq_record) < 100 and len(seq_record) % 3 == 0:
        if any(x in seq_record.description.split(" ")[3] for x in keywords):  # 3: Type
            aa_sequences += [str(seq_record.seq.translate())]
            names += [seq_record.description.split(" ")[0]]  # 0: Part name
            descriptions += [seq_record.description]


# Prepare reference dataframe
data_dict = {
    "name": pandas.Series(names),
    "sequence": pandas.Series(aa_sequences),
    "description": pandas.Series(descriptions),
}
igem_dataset = pandas.DataFrame(data=data_dict)
igem_dataset.drop_duplicates(subset=["sequence"], inplace=True, ignore_index=True)

print(len(igem_dataset))
# igem_dataset.to_csv("igem.csv", index=False)

# Alternatively, the seq_record.description can be filtered for "tag" or "linker" etc
for index, row in igem_dataset.iterrows():
    if "Tag" in row["description"]:
        print(row["description"])

# Get all part types:
part_types = []
for seq_record in SeqIO.parse("All_Parts.fasta", "fasta"):
    part_type = seq_record.description.split(" ")[3]
    part_types += [part_type]
part_types = set(part_types)
print(part_types)
