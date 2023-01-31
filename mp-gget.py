import gget
import pandas as pd
from Bio import SeqIO
import re

# load fasta file
records = list(SeqIO.parse("./data/contigs_over40kb.fa", "fasta"))


def find_proteins(dna_sequence, min_length=30):
    # regex for finding proteins with ATG as start codon and TAA, TAG or TGA as stop codon excluding the stop codon
    protein_regex = "(ATG)([ATGC]{3}){1,}?(TAA|TAG|TGA)"

    matches = re.finditer(protein_regex, dna_sequence)
    proteins = []
    proteins_index = []
    for match in matches:
        span = match.span()
        # exclude proteins shorter than 100 amino acids
        if span[1] - span[0] > min_length * 3:
            proteins.append(match.group()[:-3])
            proteins_index.append(match.span())
    return proteins, proteins_index

def reverse_complement(dna_sequence):
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return "".join([complement[base] for base in dna_sequence[::-1]])

def translate_dna_to_protein(dna_sequence):
    codon_table = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }
    protein = ""
    for i in range(0, len(dna_sequence), 3):
        codon = dna_sequence[i : i + 3]
        if codon in codon_table:
            protein += codon_table[codon]
        else:
            protein += "X"
    return protein

def run_blast(sequence, seq_id):
    print("Running blast for sequence {}".format(sequence))
    result_handle = gget.blast(sequence, "blastp", "nr")
    print(f"Finished blast for sequence {sequence} with result {result_handle}")
    # if the result is a dataframe, save it to a file
    if isinstance(result_handle, pd.DataFrame):
        result_handle.to_csv(f"results/{seq_id}.txt", index=False)
        print(f"Saved result to file results/{seq_id}.txt")
    else:
        print(f"No results found for sequence {seq_id}")

if __name__ == '__main__':
    for record in records[0:1]:
        proteins, proteins_index = find_proteins(str(record.seq), min_length=100)
        proteins = [translate_dna_to_protein(protein) for protein in proteins]
        proteins_ids = [f"{record.id}_{i[0]}-forward" for i in proteins_index]

        reverse_complement_proteins, reverse_complement_proteins_index = find_proteins(reverse_complement(str(record.seq)), min_length=100)
        reverse_complement_proteins = [translate_dna_to_protein(protein) for protein in reverse_complement_proteins]
        reverse_complement_proteins_ids = [f"{record.id}_{i[0]}-reverse" for i in reverse_complement_proteins_index]

        proteins.extend(reverse_complement_proteins)
        proteins_ids.extend(reverse_complement_proteins_ids)
        print(proteins[-5:], proteins_ids[-5:])
        print(len(proteins), len(proteins_ids))

        for protein, protein_id in zip(proteins, proteins_ids):
            run_blast(protein, protein_id)

