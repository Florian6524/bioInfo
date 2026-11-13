import itertools
from Bio import SeqIO

def load_local_fasta(file_path):
    record = SeqIO.read(file_path, "fasta")
    return str(record.seq).upper()[:200000]  #shorten the file to the first 200k since it was too long

def reverse_complement(seq):
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]

def generate_brute_force_kmers(k=6):
    return [''.join(p) for p in itertools.product('ACGT', repeat=k)]

def search_transposons(dna_seq, kmers, max_length=150):
    transposons = []
    k = len(kmers[0])
    
    for kmer in kmers:
        rev_comp = reverse_complement(kmer)
        
        start_indices = []
        idx = dna_seq.find(kmer)
        while idx != -1:
            start_indices.append(idx)
            idx = dna_seq.find(kmer, idx + 1)
            
        if not start_indices:
            continue
            
        for start in start_indices:
            search_start = start + k
            search_end = min(start + max_length, len(dna_seq))
            
            rc_idx = dna_seq.find(rev_comp, search_start, search_end)
            while rc_idx != -1:
                end = rc_idx + k
                length = end - start
                
                transposons.append({
                    "start": start,
                    "length": length,
                    "kmer": kmer
                })
                
                rc_idx = dna_seq.find(rev_comp, rc_idx + 1, search_end)
                
    unique_transposons = []
    seen = set()
    for t in transposons:
        identifier = (t["start"], t["length"], t["kmer"])
        if identifier not in seen:
            seen.add(identifier)
            unique_transposons.append(t)
            
    unique_transposons.sort(key=lambda x: x["start"])
    return unique_transposons

if __name__ == "__main__":
    dna = load_local_fasta("sequence.fasta")
    kmers = generate_brute_force_kmers(k=6)
    detected_transposons = search_transposons(dna, kmers, max_length=150)
    
    for t in detected_transposons:
        print(f"Sequence: {t['kmer']}, Start Index: {t['start']}, Length: {t['length']}")