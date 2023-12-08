import numpy as np

def neighbor_joining(matrix, names):
    n = len(matrix)
    tree = np.zeros((2 * n - 1, 3))  # Initialize the tree matrix
    index_to_name = {i: name for i, name in enumerate(names)}
    name_to_index = {name: i for i, name in enumerate(names)}

    # Function to find the minimum element in the distance matrix
    def find_min_element(dist_matrix):
        min_val = float('inf')
        min_i, min_j = -1, -1
        for i in range(n):
            for j in range(i + 1, n):
                if dist_matrix[i, j] < min_val:
                    min_val = dist_matrix[i, j]
                    min_i, min_j = i, j
        return min_i, min_j

    # Main Neighbor-Joining algorithm
    for step in range(n - 2):
        total_dist = np.sum(matrix, axis=1)
        q_matrix = (n - 2) * matrix - total_dist[:, np.newaxis] - total_dist[np.newaxis, :]

        # Find the minimum element in the Q-matrix
        i, j = find_min_element(q_matrix)

        # Calculate the new node and update the tree matrix
        new_node = n + step
        limb_i = 0.5 * (matrix[i, j] + (total_dist[i] - total_dist[j]) / (n - 2))
        limb_j = matrix[i, j] - limb_i

        tree[new_node, 0] = i
        tree[new_node, 1] = j
        tree[new_node, 2] = limb_i

        # Update the distance matrix for the new node
        matrix = np.delete(matrix, j, axis=1)
        matrix = np.delete(matrix, j, axis=0)

        new_distances = 0.5 * (matrix[i, :] + matrix[j, :] - matrix[i, j])
        matrix = np.delete(matrix, i, axis=1)
        matrix = np.delete(matrix, i, axis=0)

        matrix = np.vstack([matrix, new_distances])
        matrix = np.column_stack([matrix, new_distances.reshape(-1, 1)])

        names.append(f"Internal_{new_node + 1}")
        index_to_name[new_node] = names[-1]
        name_to_index[names[-1]] = new_node

    # Add the final two nodes to the tree
    tree[-1, 0] = name_to_index[names[-2]]
    tree[-1, 1] = name_to_index[names[-1]]
    tree[-1, 2] = 0.5 * matrix[name_to_index[names[-2]], name_to_index[names[-1]]]

    return tree, index_to_name

def write_tree_to_file(tree, output_file, index_to_name):
    with open(output_file, "w") as file:
        for row in tree:
            ancestor, descendant, length = row
            file.write(f"{index_to_name[int(ancestor)]}\t{index_to_name[int(descendant)]}\t{length}\n")

def main():
    # Example distance matrix and sequence names
    matrix = np.array([
        [0, 5, 9, 9],
        [5, 0, 10, 10],
        [9, 10, 0, 8],
        [9, 10, 8, 0]
    ])

    names = ['Seq1', 'Seq2', 'Seq3', 'Seq4']

    # Perform Neighbor-Joining
    tree, index_to_name = neighbor_joining(matrix, names)

    # Write the tree to the output file
    output_file = "edges.txt"
    write_tree_to_file(tree, output_file, index_to_name)

    print(f"Phylogenetic tree written to {output_file}")

if __name__ == "__main__":
    main()

