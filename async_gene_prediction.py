import asyncio
import gget
from Bio import SeqIO
import re

# load fasta file
records = list(SeqIO.parse("./data/contigs_over40kb.fa", "fasta"))


def find_proteins(dna_sequence, min_length=100):
    protein_regex = "(ATG)([ATGC]{3}){" + str(min_length - 1) + ",}?(TAA|TAG|TGA)" #regex for finding proteins
    matches = re.finditer(protein_regex, dna_sequence)
    proteins = []
    proteins_index = []
    for match in matches:
        proteins.append(match.group())
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


def blast_protein(protein):
    print("Blasting protein")
    return gget.blast(protein, "blastp", "nr")


async def async_wrapper(args):
    loop = asyncio.get_running_loop()
    data = await loop.run_in_executor(None, blast_protein, args)
    return data

async def async_run(proteins, func=async_wrapper):
    # Create a list of tasks to perform asynchronously
    tasks = []
    for protein in proteins:
        tasks.append(asyncio.create_task(func(protein)))

    # Wait for all tasks to complete
    results = await asyncio.gather(*tasks)
        
    return results

def main():
    proteins, proteins_index = find_proteins(str(records[0].seq))
    proteins = [translate_dna_to_protein(protein) for protein in proteins]
    
    loop = asyncio.new_event_loop()
    blast_results = loop.run_until_complete(async_run(proteins[:10], func=async_wrapper))
    return blast_results


if __name__ == "__main__":
    main()