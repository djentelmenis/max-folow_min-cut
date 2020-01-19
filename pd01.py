from collections import defaultdict
from copy import deepcopy
from re import findall


# Edmonds-Karp algorithm
# O(V * E^2)
def get_max_flow(graph, s, t):
    n = len(graph)
    r_graph = [[0] * n for _ in range(n)]
    path = bfs(graph, r_graph, s, t)
    while path is not None:
        flow = min(graph[u][v] - r_graph[u][v] for u, v in path)
        for u, v in path:
            r_graph[u][v] += flow
            r_graph[v][u] -= flow
        path = bfs(graph, r_graph, s, t)
    return sum(r_graph[s][i] for i in range(n))


# Finds shortest path by using BFS
# O(V * E)
def bfs(graph, r_graph, s, t):
    queue = [s]
    paths = {s: []}
    if s == t:
        return paths[s]
    while queue:
        u = queue.pop(0)
        for v in range(len(graph)):
            if (graph[u][v] - r_graph[u][v] > 0) and v not in paths:
                paths[v] = paths[u] + [(u, v)]
                if v == t:
                    return paths[v]
                queue.append(v)
    return None


# Finds the minimum articulation points, that would disconnect the source from the sink, by looping through
# the graph once and removing vertices in O(V), removing and restoring nodes in O(V) and using Edmonds-Karp
# in O(V * E^2) to check if that changes the max flow, so total O(V * (V + V + V * E^2)) = O(V^2 * E^2)
# Could be optimized by removing the deepcopy from the while loop
def get_min_cut(a_matrix, a_list, k, s, t):
    cut_nodes = []
    remaining_flow = k
    v = 0
    while remaining_flow > 0 and v < len(a_matrix):
        # Skips source and sink
        if v == s or v == t:
            v += 2
            continue

        # Remove a node from the graph but saves the connections for later
        for i, ref in enumerate(a_matrix[v]):
            if ref > 0:
                a_matrix[i][v] = 0
        removed_node = deepcopy(a_matrix[v])
        a_matrix[v] = [0] * len(a_matrix)

        new_flow = get_max_flow(a_matrix, s, t)
        if new_flow < remaining_flow:
            # Get name or original node by
            # finding the removed node in adjacency list by index and removing the appended name
            cut_nodes.append(list(a_list.keys())[v][:-3])
            remaining_flow = new_flow
        else:
            # if removed node had no effect, the node and connections to it are
            # restored for next iterations to avoid false positives
            a_matrix[v] = deepcopy(removed_node)
            for i in removed_node:
                if i > 0:
                    a_matrix[i][v] = 1
        v += 2

    return cut_nodes


# Creates an adjacency matrix from an adjacency list
# by traversing each node in the adjacency list in O(V^2)
# assuming the index lookup takes constant time
def create_adj_matrix(a_list):
    a_matrix = [[0 for _ in range(len(a_list))] for _ in range(len(a_list))]
    for i, n in enumerate(a_list):
        for j, m in enumerate(a_list[n]):
            # Get index for the connection the current node has by
            # finding the index in the adjacency list by the connections name
            index = list(a_list.keys()).index(m)
            a_matrix[i][index] = 1

    return a_matrix


# Creates a transformed graph G1 from G so that max-flow will find vertex-disjoint paths
# For each node V in G, create two new nodes V_In and V_Out with an edge between them in G1
# For each edge U,V in G, create an edge U_Out,V_In, and V_Out,U_In if edge undirected
# O(V+E)
def transform_graph(a_list, s, t):
    t_list = defaultdict(list)

    # convert vertexes
    for v in a_list:
        t_list[v + '_In'] = [v + '_Ou']
        t_list[v + '_Ou'] = []

    # convert edges
    for v in a_list:
        for r in a_list[v]:
            t_list[v + '_Ou'].append(r + '_In')

    return t_list


# Read input file which is at most 6*E long
# and creates an adjacency list in O(E) time
def read_file(name):
    with open(name, 'r') as f:
        input_graph = f.readlines()[0]

    # Split input in to edges
    # Not sure how findall is implemented, but it could be replaced with a loop that
    # parses the string once in O(E*4) and creates edges with dynamic programming in constant time
    edges = findall("\d+\s\d+", input_graph)

    # Creates a graph adjacency list by looping through the edges in O(E) time
    # assuming findall takes constant time by appending connections to the adjacency list
    # findall could be replaced with a O(5) loop and dynamic programming (with E < 100)
    a_list = defaultdict(list)
    for edge in edges:
        # Split edge in to nodes
        nodes = findall("\d+", edge)
        from_node = nodes[0]
        to_node = nodes[1]
        a_list[str(from_node)].append(str(to_node))
        a_list[str(to_node)].append(str(from_node))

    return a_list


# Total complexity O(V^3 * E^2)
def main(i=2, source=0, sink=20):
    adj_list = read_file(str(i) + '.txt')                                               # O(E)
    trans_list = transform_graph(adj_list, source, sink)                                # O(E+V)
    trans_adj_matrix = create_adj_matrix(trans_list)                                    # O(V^2)

    S = T = None
    # To check if requested nodes ar part of the graph
    try:
        S = list(trans_list.keys()).index(str(source) + '_Ou')
        T = list(trans_list.keys()).index(str(sink) + '_In')
    except:
        print('Source or sink outside the graph')

    if S is not None and T is not None:
        max_flow = get_max_flow(trans_adj_matrix, S, T)                                     # O(E^2 * V)
        articulation_points = get_min_cut(trans_adj_matrix, trans_list, max_flow, S, T)     # O(E^2 * V^2)

        print('There exist', max_flow, 'disjoint paths between', source, 'and', sink, 'and')
        print('There exist', len(articulation_points), 'articulation points between', source, 'and', sink, ':')
        print(articulation_points)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Finds minimum articulation points in a graph,'
                                                 'that would disconnect the source from the sink')
    parser.add_argument('input', help='which input file to run [1-3]')
    parser.add_argument('source', help='source')
    parser.add_argument('sink', help='sink')
    args = parser.parse_args()

    # Executes script with variables
    main(i=args.input, source=args.source, sink=args.sink)
