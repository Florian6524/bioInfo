import urllib.request
import random

def get_ncbi_sequence():
    """Fetch DNA sequence from NCBI or use a fallback sequence."""
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=J02459.1&retmode=text&rettype=fasta&seq_start=1&seq_stop=2000"
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=5) as response:
            fasta_data = response.read().decode('utf-8')
            lines = fasta_data.splitlines()
            seq = "".join([line.strip() for line in lines if not line.startswith(">")])
            if len(seq) >= 1000:
                return seq.upper()
    except Exception:
        pass
    
    #Fallback: just generate a DNA sequence of about 2000 bases, in case it can't fetch
    random.seed(101)
    base_seq = "".join(random.choices(['A', 'C', 'G', 'T'], k=1900))
    return base_seq[:100] + "ATGCATGC" + base_seq[108:500] + "GATCGATC" + base_seq[508:]

def ex2():
    dna_seq = get_ncbi_sequence()
    dna_seq = dna_seq[:2000]
    print(f"Analyzing DNA sequence of length: {len(dna_seq)} bp\n")
    
    detected_repeats = []
    i = 0
    n = len(dna_seq)
    
    while i < n:
        found_match = False
        for u_len in range(3, 6):
            if i + u_len > n:
                continue
            unit = dna_seq[i:i + u_len]
            
            k = 1
            while True:
                next_start = i + (k * u_len)
                if next_start + u_len > n:
                    break
                if dna_seq[next_start:next_start + u_len] == unit:
                    k += 1
                else:
                    break
            
            total_rep_len = u_len * k
            
            if 6 <= total_rep_len <= 10:
                matched_string = unit * k
                detected_repeats.append({
                    'position': i,
                    'unit': unit,
                    'unit_length': u_len,
                    'repeats': k,
                    'total_length': total_rep_len,
                    'sequence': matched_string
                })
                i += total_rep_len - 1
                found_match = True
                break
        
        i += 1

    print(f"Found {len(detected_repeats)} matching repetition(s):")
    print(f"{'Position':<10} | {'Unit':<8} | {'Unit Len':<10} | {'Repeats':<8} | {'Total Len':<10} | {'Sequence'}")
    print("-" * 70)
    
    for rep in detected_repeats:
        print(f"{rep['position']:<10} | {rep['unit']:<8} | {rep['unit_length']:<10} | {rep['repeats']:<8} | {rep['total_length']:<10} | {rep['sequence']}")

if __name__ == "__main__":
    ex2()