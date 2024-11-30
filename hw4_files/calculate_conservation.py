
with open("Homework4-seqs-with-primers.fna", "r") as file:
    lines = file.readlines()
    alignment = lines[1::2]
    file.close()

# Calculate conservation rate for each position
conservation_rates = []
alignment_length = len(alignment[0])
reads_num = len(alignment)

for i in range(alignment_length):
    non_gap_chars = [seq[i] for seq in alignment if seq[i] in "ACGT"]
    total_chars = len(non_gap_chars)

    if total_chars == 0:
        conservation_rate = 0.0
    else:
        most_common_char_count = non_gap_chars.count(max(set(non_gap_chars), key=non_gap_chars.count))
        conservation_rate = most_common_char_count / reads_num

    conservation_rates.append(conservation_rate)

# Write conservation rates to a text file
output_file_path = "solution-problem-1.txt"
with open(output_file_path, "w") as file:
    for rate in conservation_rates:
        file.write(f"{rate:.4f}\n")
    file.close()

print(f"Conservation rates saved to {output_file_path}")

