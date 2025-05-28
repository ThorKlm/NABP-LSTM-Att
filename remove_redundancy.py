from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

def sequence_identity(seq1, seq2):
    """Fast identity score between two sequences."""
    if len(seq1) != len(seq2):
        # Trim longer to match shorter length
        min_len = min(len(seq1), len(seq2))
        seq1 = seq1[:min_len]
        seq2 = seq2[:min_len]
    matches = sum(res1 == res2 for res1, res2 in zip(seq1, seq2))
    return matches / len(seq1)

def remove_redundant_nanobodies(input_fasta="data_prev/asFasta/nanobody_seqs.fasta",
                                 output_fasta="data_prev/asFasta/nanobody_seqs_98_nocdhit.fasta",
                                 threshold=0.98):
    input_records = list(SeqIO.parse(input_fasta, "fasta"))
    unique_records = []

    for i, rec in enumerate(input_records):
        is_unique = True
        for uniq in unique_records:
            identity = sequence_identity(str(rec.seq), str(uniq.seq))
            if identity >= threshold:
                is_unique = False
                break
        if is_unique:
            unique_records.append(rec)

    SeqIO.write(unique_records, output_fasta, "fasta")
    print(f"Filtered {len(input_records)} → {len(unique_records)} non-redundant sequences (threshold: {threshold * 100:.0f}%)")
    return output_fasta

# Example usage
remove_redundant_nanobodies()
