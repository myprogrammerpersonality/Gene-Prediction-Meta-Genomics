from Bio import SeqIO
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP

# load validated genes
records = list(SeqIO.parse("./results/validated_genes.fa", "fasta"))

# calculate isoelectric point for each gene
records_pi = [IP(record.seq).pi() for record in records]

# keep genes with high negative charge
records_pi_filtered = [record for record, pi in zip(records, records_pi) if pi < 5]
print(len(records_pi_filtered))

# save filtered genes
SeqIO.write(records_pi_filtered, "./results/validated_genes_potential_AMP_targets.fa", "fasta")