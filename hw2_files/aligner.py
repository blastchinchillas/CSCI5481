import argparse

def needleman_wunsch(seq1_file, seq2_file, output, gap_penalty, mismatch_penalty, match_score, unpenalized_start_end, affine_gap_penalty):
    # Initialize the alignment matrix
    with open('{}'.format(seq1_file), 'r') as file:
        lines = file.readlines()
        anno1 = lines[0].replace("\n",'')
        seq1 = lines[1].replace("\n",'')
        file.close()
    with open('{}'.format(seq2_file), 'r') as file:
        lines = file.readlines()
        anno2 = lines[0].replace("\n",'')
        seq2 = lines[1].replace("\n",'')
        file.close()
    m, n = len(seq1), len(seq2)
    alignment_matrix = [[0] * (n + 1) for _ in range(m + 1)]
    event_matrix = [["m"] * (n + 1) for _ in range(m + 1)]

    # Initialize the matrix with gap penalties
    if affine_gap_penalty == 0:
        affine_gap_penalty = gap_penalty
    for i in range(m + 1):
        alignment_matrix[i][0] = 0 if unpenalized_start_end or i == 0 else affine_gap_penalty + (i-1)*gap_penalty
        event_matrix[i][0] = "d"
    for j in range(n + 1):
        alignment_matrix[0][j] = 0 if unpenalized_start_end or j == 0 else affine_gap_penalty + (j-1)*gap_penalty
        event_matrix[0][j] = "i"
    event_matrix[0][0] = "m"

    # Fill in the alignment matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = alignment_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
            delete = alignment_matrix[i - 1][j] + (0 if j == n and unpenalized_start_end else gap_penalty if event_matrix[i-1][j] == "d" else affine_gap_penalty)
            insert = alignment_matrix[i][j - 1] + (0 if i == m and unpenalized_start_end else gap_penalty if event_matrix[i][j-1] == "i" else affine_gap_penalty)
            alignment_matrix[i][j] = max(match, delete, insert)
            events = {'m': match, 'd': delete, 'i': insert}
            event_matrix[i][j] = max(events, key=events.get)

    # Traceback to find the aligned sequences
    aligned_seq1, aligned_seq2, aligned_visual = [], [], []
    i, j = m, n
    while i > 0 and j > 0:
        current_score = alignment_matrix[i][j]
        if event_matrix[i][j] == 'm':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            aligned_visual.append('|' if seq1[i - 1] == seq2[j - 1] else 'x')
            i -= 1
            j -= 1
        elif event_matrix[i][j] == 'd':
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('_')
            aligned_visual.append(' ')
            i -= 1
        else:
            aligned_seq1.append('_')
            aligned_seq2.append(seq2[j - 1])
            aligned_visual.append(' ')
            j -= 1

    # Finish tracing back if there are remaining elements
    while i > 0:
        aligned_seq1.append(seq1[i - 1])
        aligned_seq2.append('_')
        aligned_visual.append(' ')
        i -= 1
    while j > 0:
        aligned_seq1.append('_')
        aligned_seq2.append(seq2[j - 1])
        aligned_visual.append(' ')
        j -= 1

    # Reverse the aligned sequences
    aligned_seq1.reverse()
    aligned_seq2.reverse()
    aligned_visual.reverse()

    return ''.join(anno1), ''.join(anno2), ''.join(aligned_seq1), ''.join(aligned_seq2), ''.join(aligned_visual), alignment_matrix[m][n], ''.join(output)

def main():
    parser = argparse.ArgumentParser(description='Needleman-Wunsch Global Sequence Alignment')
    parser.add_argument('-q', '--query', required=True, type=str, help='First input sequence')
    parser.add_argument('-r', '--reference', required=True, type=str, help='Second input sequence')
    parser.add_argument('-o', '--output', required=True, type=str, help='Output file')
    parser.add_argument('-g', '--gap_penalty', type=int, default=-2, help='Gap penalty (default: -2)')
    parser.add_argument('-p', '--mismatch_penalty', type=int, default=-1, help='Mismatch penalty (default: -1)')
    parser.add_argument('-m', '--match_score', type=int, default=1, help='Match score (default: 1)')
    parser.add_argument('--ignore_outer_gaps', action='store_true', help='Unpenalized start and end for both sequences')
    parser.add_argument('-s', '--affine_gap_penalty', required=False, default=0, type=int)
    args = parser.parse_args()
    
    annotation1, annotation2, aligned_seq1, aligned_seq2, visualization, aligned_score, output_file = needleman_wunsch(args.query, args.reference, args.output, args.gap_penalty, args.mismatch_penalty, args.match_score, args.ignore_outer_gaps, args.affine_gap_penalty)

    lines_to_write = ["{}\n".format(aligned_score), "{}\n".format(annotation1), "{}\n".format(aligned_seq1), "{}\n".format(visualization), "{}\n".format(aligned_seq2), "{}\n".format(annotation2)]
    with open('{}'.format(output_file), 'w') as file:
        file.writelines(lines_to_write)
        file.close()
if __name__ == "__main__":
    main()

