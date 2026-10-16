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
    'TAC':'Yeast', 'TAT':'Tyrosine', 'TAA':'Stop', 'TAG':'Tyrosine',
    'TGC':'Cysteine', 'TGT':'Cysteine', 'TGA':'Stop', 'TGG':'Tryptophan'
}

#searched with AI a list of foods for each possible amino acid
aminoAcidFoods = {
    'Alanine': 'Meat, poultry, fish, eggs, dairy, avocados, and legumes.',
    'Arginine': 'Nuts, seeds (pumpkin, sesame), turkey, chicken, pork, dairy, and seaweed.',
    'Asparagine': 'Dairy, whey, beef, poultry, eggs, fish, asparagus, and potatoes.',
    'Aspartic Acid': 'Sprouted seeds, oats, avocado, asparagus, and sugarcane.',
    'Cysteine': 'Poultry, pork, beef, eggs, dairy, garlic, onions, and broccoli.',
    'Glutamic Acid': 'Eggs, poultry, meat, fish, dairy, mushrooms, and tomatoes.',
    'Glutamine': 'Beef, chicken, fish, dairy, eggs, cabbage, and spinach.',
    'Glycine': 'Bone broth, gelatin, skin/cartilage meat, poultry, fish, and dairy.',
    'Histidine': 'Meat, poultry, fish (tuna, salmon), dairy, nuts, and whole grains.',
    'Isoleucine': 'Meat, poultry, fish, eggs, cheese, seeds, nuts, and lentils.',
    'Leucine': 'Chicken, beef, pork, fish, milk, cheese, eggs, peanuts, and brown rice.',
    'Lysine': 'Red meat, pork, poultry, parmesan cheese, cod, sardines, and eggs.',
    'Methionine': 'Sesame seeds, Brazil nuts, fish, meat, poultry, eggs, and dairy.',
    'Phenylalanine': 'Meat, fish, poultry, eggs, cheese, soybeans, nuts, and seeds.',
    'Proline': 'Bone broth, gelatin, meat, dairy, egg whites, and asparagus.',
    'Serine': 'Soybeans, eggs, seaweed, lamb, beef, chicken, and nuts.',
    'Threonine': 'Cottage cheese, poultry, fish, meat, lentils, and sesame seeds.',
    'Tryptophan': 'Turkey, chicken, milk, cheese, yogurt, fish, eggs, and oats.',
    'Tyrosine': 'Cheese, soy products, chicken, turkey, fish, peanuts, and almonds.',
    'Valine': 'Dairy, meat, soy products, peanuts, and mushrooms.'
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
    
    labels = [f"{codon}\n({codonTable.get(codon.replace('U', 'T'), 'Unknown')})" for codon in codons]

    plt.figure(figsize=(11, 6))
    bars = plt.bar(labels, frequencies, color='blue', alpha=0.85)
    
    plt.xlabel("RNA Codon & Amino Acid", fontsize=12, fontweight='bold')
    plt.ylabel("Frequency (Count)", fontsize=12, fontweight='bold')
    plt.title("Top 10 Most Frequent RNA Codons in Genome", fontsize=14, fontweight='bold')
    plt.grid(axis='y', linestyle='--', alpha=0.7)

    for bar in bars:
        height = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2., height + (max(frequencies)*0.01),
                 f"{height:,}",
                 ha='center', va='bottom', fontsize=9, rotation=0)

    plt.tight_layout()
    plt.show()

def print_dietary_recommendations(top_10):
    """Analyzes amino acids in the top 10 codons and outputs foods to avoid."""
    unique_aas = set()
    for codon, _ in top_10:
        dna_codon = codon.replace('U', 'T')
        aa = codonTable.get(dna_codon, 'Unknown')
        if aa != 'Stop' and aa != 'Unknown':
            unique_aas.add(aa)

    print("Foods reccomended to avoid for each amino acid:\n")

    for aa in sorted(unique_aas):
        foods = aminoAcidFoods.get(aa, "General protein sources.")
        print(f" • {aa.upper()}:")
        print(f"   Avoid/Limit: {foods}\n")

def main():
    fasta_file = "GCF_000591035.1_ASM59103v1_cds_from_genomic.fna"  
    
    print(f"Reading genome from {fasta_file}...")
    genome_seq = read_fasta(fasta_file)
    
    if not genome_seq:
        return
    
    codon_counts = analyze_genome_codons(genome_seq, window_size=3)
    top_10 = sorted(codon_counts.items(), key=lambda item: item[1], reverse=True)[:10]
    
    print("\n" + "="*65)
    print(f"{'Rank':<6}{'RNA Codon':<12}{'Amino Acid':<18}{'Frequency':<10}")
    print("="*65)
    
    for rank, (codon, freq) in enumerate(top_10, start=1):
        dna_key = codon.replace('U', 'T')
        amino_acid = codonTable.get(dna_key, 'Unknown')
        print(f"{rank:<6}{codon:<12}{amino_acid:<18}{freq:<10,}")
        
    print("="*65)
    
    print_dietary_recommendations(top_10)
    
    print("\nGenerating plot for top 10 RNA codons...")
    plot_top_codons(top_10)

if __name__ == "__main__":
    main()