import pandas as pd
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
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

if __name__ == "__main__":
    proteins, proteins_index = find_proteins(str(records[0].seq), min_length=100)
    proteins = [translate_dna_to_protein(protein) for protein in proteins]
    proteins_ids = [f"{records[0].id}_{i[0]}-forward" for i in proteins_index]

    sequences = []
    sequences_headers = []
    for i, index in enumerate(proteins_index):
        # if the file exist add the protein to the file otherwise skip
        if os.path.exists(f"./results/k141_3522370_{index[0]}-forward.txt"):
            df = pd.read_csv(f"./results/k141_3522370_{index[0]}-forward.txt")
            df = df.sort_values(by="E value").reset_index(drop=True)
            min_evalue = df["E value"][0]
            min_evalue_description = df["Description"][0]
            min_evalue_sciname = df["Scientific Name"][0]
            if min_evalue < 0.01:
                header = f"{proteins_ids[i]}|{min_evalue}|{min_evalue_description}|{min_evalue_sciname}"
                sequences_headers.append(header)
                sequences.append(proteins[i])

    # save the sequences and headers to biopython records
    records = []
    for i, sequence in enumerate(sequences):
        record = SeqRecord(Seq(sequence), id=sequences_headers[i], description="")
        records.append(record)

    # save the records to a fasta file
    SeqIO.write(records, "./results/validated_genes.fa", "fasta")