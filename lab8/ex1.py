import numpy as np


def needleman_wunsch(s1: str, s2: str, match: int = 2, mismatch: int = -1, gap: int = -2):
    n, m = len(s1), len(s2)

    score_matrix = np.zeros((n + 1, m + 1), dtype=int)
    for i in range(1, n + 1):
        score_matrix[i, 0] = i * gap
    for j in range(1, m + 1):
        score_matrix[0, j] = j * gap

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            sub = match if s1[i - 1] == s2[j - 1] else mismatch
            diag = score_matrix[i - 1, j - 1] + sub
            up = score_matrix[i - 1, j] + gap
            left = score_matrix[i, j - 1] + gap
            score_matrix[i, j] = max(diag, up, left)

    aln1_rev = []
    aln2_rev = []
    directions_taken = []
    step_options = []

    i, j = n, m
    while i > 0 or j > 0:
        curr = score_matrix[i, j]
        options = []

        diag_val = (
            score_matrix[i - 1, j - 1] + (match if s1[i - 1] == s2[j - 1] else mismatch)
            if i > 0 and j > 0
            else None
        )
        up_val = score_matrix[i - 1, j] + gap if i > 0 else None
        left_val = score_matrix[i, j - 1] + gap if j > 0 else None

        if diag_val == curr:
            options.append("Diagonal")
        if up_val == curr:
            options.append("Up")
        if left_val == curr:
            options.append("Left")

        step_options.append(f"({i},{j})->{options}")

        chosen = options[0]
        directions_taken.append(chosen)

        if chosen == "Diagonal":
            aln1_rev.append(s1[i - 1])
            aln2_rev.append(s2[j - 1])
            i -= 1
            j -= 1
        elif chosen == "Up":
            aln1_rev.append(s1[i - 1])
            aln2_rev.append("_")
            i -= 1
        else:
            aln1_rev.append("_")
            aln2_rev.append(s2[j - 1])
            j -= 1

    aligned_s1 = "".join(reversed(aln1_rev))
    aligned_s2 = "".join(reversed(aln2_rev))
    return aligned_s1, aligned_s2, int(score_matrix[n, m]), score_matrix, directions_taken, step_options


def print_report(s1, s2, aln1, aln2, score, matrix, directions, step_options):
    print(f"s1: {s1}")
    print(f"s2: {s2}\n")

    print("Scoring Matrix:")
    print("      -  " + "  ".join(f"{c:>3}" for c in s2))
    rows = "-" + s1
    for idx, row in enumerate(matrix):
        print(f" {rows[idx]}  " + " ".join(f"{v:>4}" for v in row))

    #print("\nTraceback Directions (bottom-right to top-left):")
    #print(" -> ".join(directions)) #print here just the directions taken, without the cell indexes

    print("\nValid moves per cell during traceback:")
    print(" | ".join(step_options))

    match_line = "".join("|" if a == b and a != "_" else " " for a, b in zip(aln1, aln2))
    print("\nFinal Alignment:")
    print(f"seq1: {aln1}")
    print(f"      {match_line}")
    print(f"seq2: {aln2}")
    print(f"\nOptimal Score: {score}\n")


if __name__ == "__main__":
    s1 = "ACCGTGAAGCCAATAC"
    s2 = "AGCGTGCAGCCAATAC"
    res = needleman_wunsch(s1, s2, match=2, mismatch=-1, gap=-2) #problem with this sequences lets say is that at traceback we only get diagonal chosen
    print_report(s1, s2, *res)
