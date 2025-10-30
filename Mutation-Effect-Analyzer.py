#!/usr/bin/env python3

import collections

fasta_path   = "bacillus_3610.fa"
ptt_path     = "bacillus.ptt"
mpileup_path = "bacillus.mpileup"


with open(fasta_path) as f:
    f.readline()
    genome = "".join(line.strip() for line in f)

genes = []
with open(ptt_path) as f:
    for line in f:
        parts = line.split("\t")
        loc   = parts[0]
        start, stop = loc.split("..")
        genes.append((int(start), int(stop)))


mutations = []
with open(mpileup_path) as f:
    for line in f:
        fields = line.split("\t")
        pos    = int(fields[1])
        ref    = fields[2]
        reads  = fields[4]
        counts = collections.Counter(ch.upper() for ch in reads if ch.upper() in "ACGT")
        counts.pop(ref, None)
        if not counts:
            continue
        alt = counts.most_common(1)[0][0]
        mutations.append((pos, ref, alt))

codon_table = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

def translate(seq):
    """Translate a nucleotide string (first frame) into protein."""
    protein = ""
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        protein += codon_table.get(codon, 'X')
    return protein

output = []
for idx, (pos, ref, alt) in enumerate(mutations, start=1):
    for start, stop in genes:
        if start <= pos <= stop:
            gene_seq = genome[start-1:stop]
            rel_index = pos - start
            mut_seq = gene_seq[:rel_index] + alt + gene_seq[rel_index+1:]
            prot_orig = translate(gene_seq)
            prot_mut  = translate(mut_seq)
            output.extend([
                f">{idx}a_nucleotide", gene_seq,
                f">{idx}b_nucleotide", mut_seq,
                f">{idx}a_translated",  prot_orig,
                f">{idx}b_translated",  prot_mut
            ])
            break

with open("finalOutput.txt", "w") as out:
    out.write("\n".join(output))