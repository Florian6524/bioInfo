geneticCodes = {

    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UAA': 'Stop', 'UAG': 'Stop',
    'UGU': 'Cys', 'UGC': 'Cys', 'UGA': 'Stop', 'UGG': 'Trp',

    'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu', 'CUG': 'Leu',
    'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro', 'CCG': 'Pro',
    'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln', 'CAG': 'Gln',
    'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg', 'CGG': 'Arg',

    'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile', 'AUG': 'Met',
    'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr', 'ACG': 'Thr',
    'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys', 'AAG': 'Lys',
    'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg', 'AGG': 'Arg',

    'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val', 'GUG': 'Val',
    'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala', 'GCG': 'Ala',
    'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu', 'GAG': 'Glu',
    'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly', 'GGG': 'Gly'
}

def translateSequence(rnaSequence):
    rnaSequence = rnaSequence.strip().upper().replace('T', 'U')
    protein = []
    startIndex = 0
    for i in range(0, len(rnaSequence) - 2):
        if rnaSequence[i:i+3] == 'AUG':
            startIndex = i
            break
    print("Start index:", startIndex)
    for i in range(startIndex, len(rnaSequence), 3):
        code = rnaSequence[i:i+3]
        if len(code) == 3:
            aminoAcid = geneticCodes.get(code, '')
            protein.append(aminoAcid) #just to print stop also, can be moved after the if-break
            if aminoAcid == 'Stop':
                break
    return '-'.join(protein)


rnaSequence = input("Enter RNA sequence: ")
proteinSequence = translateSequence(rnaSequence)
print("Protein sequence:", proteinSequence)

#input: tagcAtgaTTAgCGGcGTtTgATattaC
#expected output: Met-Ile-Ser-Gly-Val-Stop