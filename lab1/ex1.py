S = "ATCGCGTCGATCGATGTACG"
fA = 0
fT= 0
fC = 0
fG = 0

for i in range(0, len(S)):
    if S[i] == "A":
        fA += 1
    elif S[i] == "T":
        fT += 1
    elif S[i] == "C":
        fC += 1
    elif S[i] == "G":
        fG += 1

print("Frequency of A:", fA)
print("Frequency of T:", fT)  
print("Frequency of C:", fC)
print("Frequency of G:", fG)
