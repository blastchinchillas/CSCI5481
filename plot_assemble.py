from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import sys

def align_sequences(sequence_a, sequence_b):
    alignments = pairwise2.align.globalxx(sequence_a, sequence_b, one_alignment_only=True, score_only=True)
    alignment = pairwise2.align.globalxx(sequence_a, sequence_b, one_alignment_only=True)[0]
    #alignment_str = format_alignment(*alignment)
    print(format_alignment(*alignment))
    #with open('alignment_result.txt', 'w') as result_file:
        #result_file.write(alignment_str)
        #result_file.close()
    return alignment

def plot_alignment_coordinates(sequence_a, sequence_b, alignment):
    a_coordinates = []
    b_coordinates = []

    a_index = 0
    b_index = 0

    for a_base, b_base in zip(alignment.seqA, alignment.seqB):
        if a_base != '-':
            a_coordinates.append(a_index)
            a_index += 1
        else:
            a_coordinates.append(None)

        if b_base != '-':
            b_coordinates.append(b_index)
            b_index += 1
        else:
            b_coordinates.append(None)

    # Plotting with lines connecting the points
    #plt.scatter(a_coordinates, b_coordinates, marker='o')
    plt.plot(a_coordinates, b_coordinates, marker='o', markersize=1, linestyle='-', color='b', linewidth=0.2)
    plt.xlabel('Known Location Coordinates')
    plt.ylabel('Assembled Location Coordinates')
    plt.title('Sequence Alignment Coordinates')
    plt.xlim(left=0)  # Set x-axis minimum to 0
    plt.ylim(bottom=0)  # Set y-axis minimum to 0
    #plt.savefig('alignment_plot.png')
    plt.show()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python plot_assemble.py sequences_file; sequences_file should only contains 2 lines, each line is a pure DNA sequence.")
        sys.exit(1)
    with open(sys.argv[1], "r") as file:
        lines = file.readlines()
        file.close()
    sequence_a = lines[1]
    sequence_b = lines[0]
    alignment = align_sequences(sequence_a, sequence_b)
    plot_alignment_coordinates(sequence_a, sequence_b, alignment)
