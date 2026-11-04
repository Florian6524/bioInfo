import random
import urllib.request
import matplotlib.pyplot as plt


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
    except Exception as e:
        print(f"Notice: Could not fetch live NCBI sequence ({e}). Using generated sequence.")
    
    #Fallback: just generate a DNA sequence of about 2000 bases, in case it can't fetch
    random.seed(42)
    return "".join(random.choices(['A', 'C', 'G', 'T'], weights=[0.29, 0.21, 0.21, 0.29], k=2000))

def ex1():
    dna_seq = get_ncbi_sequence()
    seq_len = len(dna_seq)
    print(f"Retrieved DNA sequence of length: {seq_len} bp\n")
    
    samples = []
    for i in range(10):
        sample_len = random.randint(100, min(1000, seq_len - 50))
        start_pos = random.randint(0, seq_len - sample_len)
        sample = dna_seq[start_pos:start_pos + sample_len]
        samples.append((f"Sample_{i+1}", sample, len(sample)))
    
    print(f"{'Sample Name':<12} | {'Length (bp)':<12} | {'Start Pos'}")
    print("-" * 42)
    for name, seq, length in samples:
        print(f"{name:<12} | {length:<12}")
    
    print("\nSimulating Gel Electrophoresis Migration")
    print("Note: Shorter fragments experience less friction and travel further.\n")
    
    sorted_samples = sorted(samples, key=lambda x: x[2], reverse=True)
    
    max_len = max(s[2] for s in samples)
    min_len = min(s[2] for s in samples)
    
    print("(Origin / Wells)")
    print("-" * 50)
    for name, seq, length in sorted_samples:
        if max_len == min_len:
            pos = 20
        else:
            pos = int(2 + 43 * (max_len - length) / (max_len - min_len))
        lane = " " * pos + " " * max(2, int(length / 50)) + f"  ({name}: {length} bp)"
        print(lane)
    print("-" * 50)
    print("(Anode / Positive Electrode)")

    fig, ax = plt.subplots(figsize=(6, 8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, max_len * 1.1)
    ax.invert_yaxis()
        
    for idx, (name, seq, length) in enumerate(samples):
        migration_distance = max_len - length + (max_len * 0.05)
        ax.barh(migration_distance, 0.6, left=idx+0.7, height=max_len*0.03, color='indigo')
        ax.text(idx + 1, migration_distance - (max_len * 0.02), f"{length}bp", ha='center', fontsize=8)
            
    ax.set_xticks(range(1, len(samples) + 1))
    ax.set_xticklabels([s[0] for s in samples], rotation=45, ha='right')
    ax.set_ylabel("Migration Distance (bp equivalent)")
    ax.set_title("Gel Electrophoresis Simulation")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    ex1()