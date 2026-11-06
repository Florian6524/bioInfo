def find_transposons(dna_sequence, ir_length=4, min_length=10, max_length=60):
    detected_elements = []
    n = len(dna_sequence)
    
    i = 0
    while i < n - ir_length:
        left_ir = dna_sequence[i : i + ir_length]
        
        j = i + ir_length + min_length
        while j <= min(n - ir_length, i + max_length):
            right_ir = dna_sequence[j : j + ir_length]
            
            if right_ir == left_ir:
                start_pos = i
                end_pos = j + ir_length
                
                detected_elements.append({
                    "start": start_pos,
                    "end": end_pos,
                    "left_ir": left_ir,
                    "right_ir": right_ir,
                    "length": end_pos - start_pos,
                    "sequence": dna_sequence[start_pos:end_pos]
                })
                break
            j += 1
        i += 1
        
    return detected_elements

if __name__ == "__main__":
    try:
        with open("ex1_dna.txt", "r") as f:
            dna_sequence = f.read().strip()
    except FileNotFoundError:
        print("Error: 'ex1_dna.txt' not found. Please run Exercise 1 first")
        exit()

    print(f"Loaded DNA Sequence (Length: {len(dna_sequence)} bases)[cite: 1]\n")
    transposons = find_transposons(dna_sequence, ir_length=4)
    
    print(f"Total Transposable Elements Detected: {len(transposons)}")
    for idx, te in enumerate(transposons, 1):
        print(f"\nTransposon {idx}")
        print(f"Start Position : Index {te['start']}")
        print(f"End Position   : Index {te['end']}")
        print(f"Length         : {te['length']} bases")
        print(f"Inverted Repeats: {te['left_ir']} ... {te['right_ir']}")
        print(f"Full Sequence  : {te['sequence']}")