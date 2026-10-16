import matplotlib.pyplot as plt

codonTable = {
    'ATA':'Isoleucine', 'ATC':'Isoleucine', 'ATT':'Isoleucine', 'ATG':'Methionine',
    'ACA':'Threonine', 'ACC':'Threonine', 'ACG':'Threonine', 'ACT':'Threonine',
    'AAC':'Asparagine', 'AAT':'Asparagine', 'AAA':'Lysine', 'AAG':'Lysine',
    'AGC':'Serine', 'AGT':'Serine', 'AGA':'Arginine', 'AGG':'Arginine',
    'CTA':'Leucine', 'CTC':'Leucine', 'CTG':'Leucine', 'CTT':'Leucine',
    'CCA':'Proline', 'CCC':'Proline', 'CCG':'Proline', 'CCT':'Proline',
    'CAC':'Histidine', 'CAT':'Histidine', 'CAA':'Glutamine', 'CAG':'Glutamine',
    'CGA':'Arginine', 'CGC':'Arginine', 'CGG':'Arginine', 'CGT':'Arginine',
    'GTA':'Valine', 'GTC':'Valine', 'GTG':'Valine', 'GTT':'Valine',
    'GCA':'Alanine', 'GCC':'Alanine', 'GCG':'Alanine', 'GCT':'Alanine',
    'GAC':'Aspartic Acid', 'GAT':'Aspartic Acid', 'GAA':'Glutamic Acid', 'GAG':'Glutamic Acid',
    'GGA':'Glycine', 'GGC':'Glycine', 'GGG':'Glycine', 'GGT':'Glycine',
    'TCA':'Serine', 'TCC':'Serine', 'TCG':'Serine', 'TCT':'Serine',
    'TTC':'Phenylalanine', 'TTT':'Phenylalanine', 'TTA':'Leucine', 'TTG':'Leucine',
    'TAC':'Tyrosine', 'TAT':'Tyrosine', 'TAA':'Stop', 'TAG':'Stop',
    'TGC':'Cysteine', 'TGT':'Cysteine', 'TGA':'Stop', 'TGG':'Tryptophan'
}

def read_fasta(file_path):
    sequence_parts = []
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('>'):
                    continue
                sequence_parts.append(line.upper())
        return "".join(sequence_parts)
    except FileNotFoundError:
        print(f"Error: The file '{file_path}' was not found.")
        return None

def analyze_genome_codons(sequence, window_size=3):
    valid_nuc = set('ATCG')
    clean_seq = "".join([n for n in sequence if n in valid_nuc])
    
    three = {}
    for i in range(len(clean_seq) - window_size + 1):
        dna_codon = clean_seq[i:i+window_size]
        rna_codon = dna_codon.replace('T', 'U')
        
        three[rna_codon] = three.get(rna_codon, 0) + 1
        
    return three

def plot_top_codons(top_10):
    codons = [item[0] for item in top_10]
    frequencies = [item[1] for item in top_10]
    
    labels = [f"{codon}\n({codonTable.get(codon, '?')})" for codon in codons]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, frequencies, color='blue', alpha=0.85)
    
    plt.xlabel("Codon (Amino Acid)", fontsize=12, fontweight='bold')
    plt.ylabel("Frequency (Count)", fontsize=12, fontweight='bold')
    plt.title("Top 10 Most Frequent Codons in Genome", fontsize=14, fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + (max(frequencies)*0.01),
                 f"{height:,}",
                 ha='center', va='bottom', fontsize=10, rotation=0)

    plt.tight_layout()
    plt.show()

def main():
    fasta_file = "GCF_000591035.1_ASM59103v1_cds_from_genomic.fna"  
    
    genome_seq = read_fasta(fasta_file)
    
    if not genome_seq:
        return
    
    codon_counts = analyze_genome_codons(genome_seq, window_size=3)
    
    top_10 = sorted(codon_counts.items(), key=lambda item: item[1], reverse=True)[:10]
    
    print("\n" + "="*50)
    print(f"{'Rank':<6}{'Codon':<8}{'Amino Acid':<14}{'Frequency':<10}")
    print("="*50)
    
    for rank, (codon, freq) in enumerate(top_10, start=1):
        amino_acid = codonTable.get(codon, 'Unknown')

        rna_codon = codon.replace('T', 'U')
        
        print(f"{rank:<6}{codon} ({rna_codon})".ljust(14) + f"{amino_acid:<14}{freq:<10,}")
        
    print("="*50)
    plot_top_codons(top_10)

if __name__ == "__main__":
    main()