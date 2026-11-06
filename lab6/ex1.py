import random

def generate_artificial_dna():
    bases = ['A', 'C', 'G', 'T']
    
    te1 = "TACG" + "AAACCGGTTGACCAGGTAAACC" + "TACG"
    te2 = "GATC" + "TTAAAGGCCTTTTAA" + "GATC"
    te3 = "CGAT" + "CCGGTAACCGGTTGACCG" + "CGAT"
    
    flank_start = "".join(random.choice(bases) for _ in range(70))
    spacer = "".join(random.choice(bases) for _ in range(50))
    flank_end = "".join(random.choice(bases) for _ in range(120))

    dna_sequence = flank_start + te1 + spacer + te2 + te3 + flank_end
    
    print(f"Generated DNA Sequence (Length: {len(dna_sequence)} bp)[cite: 1]:")
    print(dna_sequence)
    
    with open("ex1_dna.txt", "w") as f:
        f.write(dna_sequence)
        
    return dna_sequence

if __name__ == "__main__":
    generate_artificial_dna()