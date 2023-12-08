import networkx as nx
import matplotlib.pyplot as plt
import sys
import time

# We need large recursion depth for whole genome assemble
sys.setrecursionlimit(10000)

def create_de_bruijn_graph(sequences, k):
    kmers = []
    for sequence in sequences:
        kmers = kmers + [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    
    # Remove duplicated kmers
    # unique_kmers = set(kmers)
    # kmers = list(unique_kmers)
    # print(len(kmers))
    graph = nx.DiGraph()
    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        # We use "count" to mark the edge for DFS, "dup" to calculate the duplication of kmers
        if graph.has_edge(prefix, suffix):
            graph[prefix][suffix]['dup'] += 1
        else:
            graph.add_edge(prefix, suffix, count=1, dup=1)
    return graph

def find_start_node(graph):
    start_nodes = [node for node in graph.nodes() if graph.in_degree(node) == 0 and graph.out_degree(node) > 0]
    # If we can not find the start node, then we arbitraryly set one of the nodes which has the max out degree as start node
    if len(start_nodes) == 0:
        start_nodes = [max(graph.nodes(), key=graph.out_degree)]
    return start_nodes

def find_eulerian_path(graph):
    # Using DFS to find eulerian path
    def dfs(node, circuit=[]):
        nonlocal circuit_long
        circuit = circuit + [node]
        if len(circuit) > len(circuit_long):
            circuit_long = circuit
        successors = list(graph.successors(node))
        for successor in successors:
            if graph[node][successor]['count'] == 1:
                graph[node][successor]['count'] = 0
                dfs(successor, circuit.copy())

    # Perform DFS to find the Eulerian path
    max_circuit = []
    start_nodes = find_start_node(graph)
    for start_node in start_nodes:
        circuit_long = []
        dfs(start_node)
        if len(circuit_long) > len(max_circuit):
            max_circuit = circuit_long
        edges = graph.edges()
        for edge in edges:
            source, target = edge
            graph[source][target]['count'] = 1
    # Assemble the genome
    assembled_genome = max_circuit[0] + ''.join([s[-1] for s in max_circuit[1:]])
    return assembled_genome

def find_weight_path(graph):
    # Find the weighted path using Kmer duplications
    def add_successor(node):
        nonlocal circuit
        circuit = circuit + [node]
        successors = list(graph.successors(node))
        valid_nodes = [x for x in successors if graph[node][x]['count'] == 1]
        if len(valid_nodes) > 0:
            successor = max(valid_nodes, key=lambda x: graph[node][x]['dup'])
            graph[node][successor]['count'] = 0
            add_successor(successor)
    
    start_nodes = find_start_node(graph)
    circuit_long = []
    for start_node in start_nodes:
        circuit = []
        add_successor(start_node)
        if len(circuit) > len(circuit_long):
            circuit_long = circuit
        edges = graph.edges()
        for edge in edges:
            source, target = edge
            graph[source][target]['count'] = 1
    # Assemble the genome
    assembled_genome = circuit_long[0] + ''.join([s[-1] for s in circuit_long[1:]])
    return assembled_genome

def visualize_graph(graph):
    pos = nx.spring_layout(graph)
    nx.draw(graph, pos, with_labels=True, font_weight='bold', node_size=700, node_color='skyblue', font_size=8)
    plt.show()

# Example usage
#sequences = ["AAGATTCTCTAC", "TCTACGA"]

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python genome_assemble.py <input fastq file> <kmer length>. \nThe de bruijn graph will be saved as 'graph.graphml'.\nTo directly visualize the graph, modify the sentence '#visualize_graph(de_bruijn_graph)'.\nTo use weight algorithm instead of DFS, modify the sentence '#assembled_contig = find_weight_path(component_graph)'.")
        sys.exit(1)
    fastq_file = sys.argv[1]
    k_value = int(sys.argv[2])

    with open(fastq_file,"r") as file:
        lines = file.readlines()
        sequences = lines[1::4]
        file.close()

    start_time1 = time.time()
    de_bruijn_graph = create_de_bruijn_graph(sequences, k_value)
    nx.write_graphml(de_bruijn_graph, "graph.graphml")
    end_time1 = time.time()
    elapsed_time1 = end_time1 - start_time1
    print(f"{elapsed_time1} second")

    start_time2 = time.time()
    assembled_genomes = []
    weakly_connected_components = list(nx.weakly_connected_components(de_bruijn_graph))
    for i in range(len(weakly_connected_components)):
        component_graph = de_bruijn_graph.subgraph(weakly_connected_components[i])
        assembled_contig = find_eulerian_path(component_graph)
        #assembled_contig = find_weight_path(component_graph)
        assembled_genomes = assembled_genomes + [assembled_contig]
    print(f"Assembled Genome:")
    sorted_list = sorted(assembled_genomes, key=len, reverse=True)
    for s in sorted_list:
        print(s)
    end_time2 = time.time()
    elapsed_time2 = end_time2 - start_time2
    print(f"{elapsed_time2} second")

    #start_time3 = time.time()
    #visualize_graph(de_bruijn_graph)
    #end_time3 = time.time()
    #elapsed_time3 = end_time3 - start_time3
    #print(f"{elapsed_time3} second")
