import matplotlib.pyplot as plt

S = "ATCGCGTCGATCGATGTACG"

two = {}
for i in range(len(S) - 1):
    x = S[i:i+2]
    two[x] = two.get(x, 0) + 1

three = {}
for i in range(len(S) - 2):
    x = S[i:i+3]
    three[x] = three.get(x, 0) + 1

plt.bar(two.keys(), two.values())
plt.title("2-letter combinations")
plt.xlabel("Combination")
plt.ylabel("Frequency")
plt.show()

plt.bar(three.keys(), three.values())
plt.title("3-letter combinations")
plt.xlabel("Combination")
plt.ylabel("Frequency")
plt.xticks(rotation=90)
plt.show()
