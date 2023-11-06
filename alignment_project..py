from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment


file1 = "COL1A1_1277.fasta"
file2 = "COL2A1_1280.fasta"

sequences_file1 = SeqIO.to_dict(SeqIO.parse(file1, "fasta"))
sequences_file2 = SeqIO.to_dict(SeqIO.parse(file2, "fasta"))

start_position = 10
end_position = 30

for id, seq1 in sequences_file1.items():
    fragment1 = seq1.seq[start_position-1:end_position]  # Adjust for 0-based indexing
    print(f"Fragment from sequence 1 ({id}): {fragment1}")

for id, seq2 in sequences_file2.items():
    fragment2 = seq2.seq[start_position-1:end_position]  # Adjust for 0-based indexing
    print(f"Fragment from sequence 2 ({id}): {fragment2}")

# the output shows that for sequence 1 (with the ID NC_000017.11:c50201631-50184101), the fragment from positions 10 to 30 is AGTTTCTCCTCGGGGTCGGAG.
# Similarly, for sequence 2 (with the ID NC_000012.12:c48006212-47972967), the fragment from positions 10 to 30 is ATAGGACCTTCTGAGCCCCAA.

# Define a function to perform alignment
def perform_alignment(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = max(alignments, key=lambda x: x[2])
    return best_alignment

# output - BiopythonDeprecationWarning: Bio.pairwise2 has been deprecated, and we intend to remove it in a future release of Biopython. As an alternative, please consider using Bio.Align.PairwiseAligner as a replacement, and contact the Biopython developers if you still need the Bio.pairwise2 module.
#   warnings.warn(
# Fragment from sequence 1 (NC_000017.11:c50201631-50184101): AGTTTCTCCTCGGGGTCGGAG
# Fragment from sequence 2 (NC_000012.12:c48006212-47972967): ATAGGACCTTCTGAGCCCCAA 

# Assuming you have already extracted the fragments
# For example:

# Extract fragments
fragment1 = "AGTTTCTCCTCGGGGTCGGAG"
fragment2 = "ATAGGACCTTCTGAGCCCCAA"

# Perform alignment
alignment_result = perform_alignment(fragment1, fragment2)
print(format_alignment(*alignment_result))

alignment_score = alignment_result[2]
identity_percentage = (alignment_result[2] / len(fragment1)) * 100

print(f"Alignment Score: {alignment_score}")
print(f"Identity Percentage: {identity_percentage}%")

# Define your sequences (replace these with your actual sequences)
sequences = [
    "AGTTTCTCCTCGGGGTCGGAG",
    "ATAGGACCTTCTGAGCCCCAA",
    # Add more sequences as needed
]
# Write the sequences to a temporary FASTA file
fasta_file = "temp.fasta"
with open(fasta_file, "w") as f:
    for i, seq in enumerate(sequences):
        f.write(f">Seq{i}\n{seq}\n")
       

#output- 

# Fragment from sequence 1 (NC_000017.11:c50201631-50184101): AGTTTCTCCTCGGGGTCGGAG
# Fragment from sequence 2 (NC_000012.12:c48006212-47972967): ATAGGACCTTCTGAGCCCCAA
# AGT------TTCT---CCTCGGGGTCGGAG-
# | |      ||||   || |     |  |  
# A-TAGGACCTTCTGAGCC-C-----C--A-A
#   Score=11
# 
# Alignment Score: 11.0
# Identity Percentage: 52.38095238095239%

# The alignment result visually shows the alignment of the two sequences with gaps represented by dashes (-). You can see where the sequences match and where they differ.
# The alignment score is 11.0, indicating the score of this alignment. Higher scores typically indicate better alignments.
# The identity percentage is approximately 52.38%, which means that about 52.38% of the positions in the alignment are identical.


# successfully implemented pairwise sequence alignment using Biopython.

