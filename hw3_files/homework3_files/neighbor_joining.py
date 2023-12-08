import sys
import numpy as np

def calculate_genetic_distance(seq1, seq2):
    length = min(len(seq1), len(seq2))
    mismatches = sum(a != b for a, b in zip(seq1, seq2))
    distance = (mismatches / length)
    return distance

def write_genetic_distances(sequences, output_file):
    with open(output_file, 'w') as f:
        # Write header
        f.write('\t'.join(sequences.keys()) + '\n')
        # Write rows
        for seq_id1, seq1 in sequences.items():
            f.write(seq_id1 + '\t')
            distance = []
            for seq_id2, seq2 in sequences.items():
                distance.append(calculate_genetic_distance(seq1, seq2))
            line = '\t'.join(map(str, distance)) + '\n'
            f.write(line)

def neighbor_joining(matrix, names):
    n = len(matrix)
    tree = np.zeros((2 * n - 3, 3))  # Initialize the tree matrix
    par = np.array(['str' for i in range(2*n -3)], dtype=object)
    # Function to find the minimum element in the distance matrix
    nodes = [i for i in range(1,len(names)+1)] # Initialize the nodes number
    def find_min_element(dist_matrix):
        min_val = float('inf')
        min_i, min_j = -1, -1
        for i in range(len(dist_matrix)):
            for j in range(i + 1, len(dist_matrix)):
                if dist_matrix[i, j] < min_val:
                    min_val = dist_matrix[i, j]
                    min_i, min_j = i, j
        return min_i, min_j
    # Main Neighbor-Joining algorithm
    for step in range(n - 2):
        total_dist = np.sum(matrix, axis=1)
        q_matrix = (n - 2 - step) * matrix - total_dist[:, np.newaxis] - total_dist[np.newaxis, :]
        # Find the minimum element in the Q-matrix
        i, j = find_min_element(q_matrix)
        # Calculate the new node and update the tree matrix
        new_node = n + step
        limb_i = 0.5 * (matrix[i, j] + ((total_dist[i] - total_dist[j]) / (n - 2 - step)))
        limb_j = matrix[i, j] - limb_i
        # Update the distance matrix for the new node
        new_distances = 0.5 * (matrix[i, :] + matrix[j, :] - matrix[i, j])
        matrix = np.column_stack([matrix, new_distances.reshape(-1, 1)])
        matrix = np.vstack([matrix, np.append(new_distances,0)])
        matrix = np.delete(matrix, j, axis=0)
        matrix = np.delete(matrix, i, axis=0)
        matrix = np.delete(matrix, j, axis=1)
        matrix = np.delete(matrix, i, axis=1)
        nodes.append(new_node+1)
        tree[2*step, 0] = nodes[-1]
        tree[2*step, 1] = nodes[i]
        tree[2*step, 2] = limb_i
        par[2*step] = str(nodes[i]) if nodes[i] <= n else ' '.join(np.char.mod('%s', par[tree[:, 0] == nodes[i]]))
        tree[2*step+1, 0] = nodes[-1]
        tree[2*step+1, 1] = nodes[j]
        tree[2*step+1, 2] = limb_j
        par[2*step+1] = str(nodes[j]) if nodes[j] <= n else ' '.join(np.char.mod('%s', par[tree[:, 0] == nodes[j]]))
        del nodes[j]
        del nodes[i]
        names.append(f"Internal_{new_node + 1}")
    # Add the final node
    tree[-1, 0] = nodes[1]
    tree[-1, 1] = nodes[0]
    tree[-1, 2] = 0.5 * matrix[0,1]
    par[-1] = str(nodes[0]) if nodes[0] <= n else ' '.join(np.char.mod('%s', par[tree[:, 0] == nodes[0]]))
    return tree,par

def write_tree_to_file(tree, output_file):
    with open(output_file, "w") as file:
        for row in tree:
            ancestor, descendant, length = row
            file.write(f"{int(ancestor)}\t{int(descendant)}\t{length}\n")
        file.close()

if __name__ == "__main__":
    if len(sys.argv) < 1:
        print("Usage: python3 neighbor_joining.py hw3.fna")
        sys.exit(1)
    fasta_file = sys.argv[1]
    output_file = "genetic-distances.txt"
    # Read sequences from FASTA file
    sequences = {}
    with open('{}'.format(fasta_file), 'r') as file:
        lines = file.readlines()
        anno = lines[0::2]
        seq = lines[1::2]
        for i in range(len(anno)):
            sequences[anno[i].replace("\n",'').replace(">",'')] = str(seq[i].replace("\n",''))
        file.close()

    # Calculate and write genetic distances
    write_genetic_distances(sequences, output_file)
    ary = []
    with open('{}'.format(output_file), 'r') as file:
        lines = [item.strip() for item in file.readlines()]
        names = lines[0].split("\t")
        for i in range(1,len(lines)):
            ary.append(lines[i].split("\t")[1:])
        matrix = np.array(ary, dtype=float)
        file.close()

    # Perform Neighbor-Joining
    tree,par = neighbor_joining(matrix, names)
    edges_file = "edges.txt"
    write_tree_to_file(tree, edges_file)

    # Calculate the bootstrap samples
    nd = [[0,0] for i in range(len(par)//2)]
    m = 100 # set pseudo-replicate times
    tmp = "tmp.txt" # set temparary file to store genetic distance 
    for i in range(m):
        sequences = {}
        with open('{}'.format(fasta_file), 'r') as file:
            lines = file.readlines()
            anno = lines[0::2]
            seq = lines[1::2]
            for j in range(len(seq)):
                seq[j] = list(str(seq[j].replace("\n",'')))
            seq_ary = np.array(seq)
            bs_ary = np.zeros_like(seq_ary)
            for j in range(len(seq_ary[0])):
                random_column_index = np.random.randint(0, seq_ary.shape[1])
                bs_ary[:, j] = seq_ary[:, random_column_index]
            for j in range(len(anno)):
                sequences[anno[j].replace("\n",'').replace(">",'')] = np.array2string(bs_ary[j], separator='')
            file.close()
        write_genetic_distances(sequences, tmp)
        ary = []
        with open('{}'.format(tmp), 'r') as file:
            lines = [item.strip() for item in file.readlines()]
            names = lines[0].split("\t")
            for j in range(1,len(lines)):
                ary.append(lines[j].split("\t")[1:])
            matrix = np.array(ary, dtype=float)
            file.close()
        tree_bs,par_bs = neighbor_joining(matrix, names)
        for j in range(len(nd)):
            nd[j][0] = len(lines) + j
            for k in range(len(nd)):
                if (set(par_bs[2*k].split(' ')) == set(par[2*j].split(' ')) and set(par_bs[2*k + 1].split(' ')) == set(par[2*j + 1].split(' '))) or (set(par_bs[2*k + 1].split(' ')) == set(par[2*j].split(' ')) and set(par_bs[2*k].split(' ')) == set(par[2*j + 1].split(' '))):
                    nd[j][1] += 0.01
                    break
    for i in range(len(nd)):
        nd[i][1] = round(nd[i][1],2)
    with open('bootstrap.txt', 'w') as file:
        for i in nd:
            line = '\t'.join(map(str, i)) + '\n'
            file.write(line)
        file.close()
        
