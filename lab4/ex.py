import random

def load_fasta(filename):
    sequences = []
    current_seq = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append("".join(current_seq))
                    current_seq = []
            else:
                current_seq.append(line.upper())
        if current_seq:
            sequences.append("".join(current_seq))
    return "".join(sequences)

def generate_reads(genome, a, b, min_coverage=10):
    genome_len = len(genome)
    coverage = [0] * genome_len
    reads = []
    
    while min(coverage) < min_coverage:
        length = random.randint(a, b)
        start = random.randint(0, max(0, genome_len - length))
        actual_len = min(length, genome_len - start)
        read = genome[start:start+actual_len]
        
        reads.append((start, read))
        for i in range(start, start + actual_len):
            coverage[i] += 1
            
    return reads, coverage

def find_overlap(s1, s2, min_length=4):
    start = 0
    while True:
        start = s1.find(s2[:min_length], start)
        if start == -1:
            return 0
        if s2.startswith(s1[start:]):
            return len(s1) - start
        start += 1

def assemble_genome(reads_collection):
    reads = list(set([r[1] for r in reads_collection]))
    current_sequences = list(reads)
    
    while len(current_sequences) > 1:
        best_overlap = -1
        best_i, best_j = -1, -1
        found = False
        
        for i in range(len(current_sequences)):
            s1 = current_sequences[i]
            for j in range(len(current_sequences)):
                if i == j:
                    continue
                s2 = current_sequences[j]
                if s2 in s1:
                    current_sequences.pop(j)
                    found = True
                    break
            if found:
                break
        
        if found:
            continue

        for i in range(len(current_sequences)):
            for j in range(len(current_sequences)):
                if i == j:
                    continue
                s1 = current_sequences[i]
                s2 = current_sequences[j]
                overlap = find_overlap(s1, s2, min_length=4)
                if overlap > best_overlap:
                    best_overlap = overlap
                    best_i, best_j = (i, j)
                    
        if best_overlap > 0:
            s1 = current_sequences[best_i]
            s2 = current_sequences[best_j]
            merged = s1 + s2[best_overlap:]
            idx1, idx2 = sorted([best_i, best_j], reverse=True)
            current_sequences.pop(idx1)
            current_sequences.pop(idx2)
            current_sequences.append(merged)
        else:
            break
            
    return current_sequences[0] if current_sequences else ""

if __name__ == "__main__":
    filename = "GCF_000591035.1_ASM59103v1_cds_from_genomic.fna"
    
    a = 100
    b = 150
    min_overlaps = 10

    print("Loading genome from file...")
    full_genome = load_fasta(filename)
    
    genome = full_genome[0:500] #used only a small part of the actual genome, it was too big, so it was taking too long
    
    print(f"Genome length: {len(genome)}")

    print("Generating reads...")
    reads, coverage_profile = generate_reads(genome, a, b, min_coverage=min_overlaps)
    print(f"Generated {len(reads)} reads. Min coverage depth: {min(coverage_profile)}")

    print("Assembling genome...")
    assembled_genome = assemble_genome(reads)

    print(f"Original:  {genome}")
    print(f"Assembled: {assembled_genome}")
    print(f"Match:     {genome == assembled_genome}")