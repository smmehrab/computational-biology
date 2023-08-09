"""
-------------------------------------------------
Submitted By:
-------------------------------------------------
Full Name   :   S.M.Mehrabul Islam
Roll        :   SH-86

-------------------------------------------------
Lab Assignment 1 [De Bruijn Fragment Assembly]
-------------------------------------------------

Given, 
Sequence = GACTTACGTACT 
k = 3

For the given input (Sequence & k),
(1) Generated the Fragments from the Sequence
        ['GAC', 'ACT', 'CTT', 'TTA', 'TAC', 'ACG', 'CGT', 'GTA', 'TAC', 'ACT']
(2) Built a De Bruijn Graph where the graph is represented by a python dictionary.
    Number of Nodes: 7
        GA: ['AC']
        AC: ['CT', 'CG']
        CT: ['TT']
        TT: ['TA']
        TA: ['AC']
        CG: ['GT']
        GT: ['TA']
(3) Check if the graph is Eulerian or not.
        Eulerian Check: [True/False]
(4) Hierholzer's Algorithm
    Finds a single eulerian path in O(V+E) time.
        ['GA', 'AC', 'CT', 'TT', 'TA', 'AC', 'CG', 'GT', 'TA', 'AC', 'CT']
    If there's only one Eulerian Path possible in the graph, 
    then this algorithm will efficiently find the path.
(5) Backtracking Algorithm
    Finds all eulerian paths via backtracking to explore all possible options.
        ['GA', 'AC', 'CT', 'TT', 'TA', 'AC', 'CG', 'GT', 'TA', 'AC', 'CT']
        ['GA', 'AC', 'CG', 'GT', 'TA', 'AC', 'CT', 'TT', 'TA', 'AC', 'CT']
(6) Got Assembled Sequences from Eulerian Paths
        "GACGTACTTACT"
        "GACTTACGTACT"

-------------------------------------------------
"""

import copy
import pydot

LINE = "-------------------------------------------------"

class FragmentGenerator:
    def __init__(self):
        pass

    def generate_fragments(self, sequence, k) :
        """
            Generate fragments from the given sequence.

            Args:
                sequence (str): input DNA sequence.
                k (int): The size of the fragments.

            Returns:
                List[str]: A list of generated fragments.
        """
        fragments = []
        sequence_length = len(sequence)
        for i in range(0, sequence_length - k + 1):
            fragment = sequence[i: i + k]
            fragments.append(fragment)
        return fragments

class DeBruijnAssembler:
    def __init__(self):
        pass

    def _find_degrees(self, graph):
        """
            Find the in-degree and out-degree of each node in the graph.

            Args:
                graph (Dict[str, List[str]]): The De Bruijn graph.

            Returns:
                Tuple[Dict[str, int], Dict[str, int]]: tuple containing dictionaries of in-degree and out-degree for each node.
        """
        indegree = {node: 0 for node in graph}
        outdegree = {node: 0 for node in graph}

        for node, neighbors in graph.items():
            for neighbor in neighbors:
                indegree[neighbor] += 1
            outdegree[node] = len(neighbors)

        return indegree, outdegree

    def _eulerian_check(self, graph):
        """
            Check if the given De Bruijn graph has an Eulerian path.

            Args:
                graph (Dict[str, List[str]]): De Bruijn graph.

            Returns:
                Tuple[bool, Optional[str]]: tuple containing a boolean indicating if the graph has an Eulerian path and the starting node for the path (if exists).
        """
        indegree, outdegree = self._find_degrees(graph)

        odd_indegree_count = 0
        odd_outdegree_count = 0
        start_node = None

        for node, neighbors in graph.items():
            if outdegree[node] - indegree[node] > 1 or indegree[node] - outdegree[node] > 1:
                return False, None
            elif outdegree[node] - indegree[node] == 1:
                odd_outdegree_count += 1
                start_node = node
            elif indegree[node] - outdegree[node] == 1:
                odd_indegree_count += 1
            elif indegree[node] != outdegree[node]:
                return False, None

        if odd_indegree_count == 1 and odd_outdegree_count == 1:
            return True, start_node
        elif odd_indegree_count == 0 and odd_outdegree_count == 0:
            for node in graph:
                if outdegree[node] > 0:
                    start_node = node
                    break
            return True, start_node

        return False, None

    def _hierholzer(self, graph, start_node):
        """
            Find an Eulerian path in the given De Bruijn graph using the Hierholzer Algorithm.

            Args:
                graph (Dict[str, List[str]]): De Bruijn graph.
                start_node (str): Starting node for the Eulerian path.

            Returns:
                List[str]: eulerian path
        """
        if len(graph) == 0:
            return []
        euler = []
        stack = [start_node]
        while stack:
            node = stack[-1]
            if graph[node]:
                next_node = graph[node].pop()
                stack.append(next_node)
            else:
                euler.append(stack.pop())
        return euler[::-1]

    def _init_edge_count(self, graph):
        """
            Initialize the edge status dictionary for the given De Bruijn graph.

            Args:
                graph (Dict[str, List[str]]): De Bruijn graph.

            Returns:
                Dict[Tuple[str, str], int]: Initialized edge status dictionary.
        """
        edge_count = {}
        for node, neighbors in graph.items():
            for neighbor in neighbors:
                if (node, neighbor) not in edge_count:
                    edge_count[(node, neighbor)] = 0
                edge_count[(node, neighbor)] += 1
        return edge_count

    def _is_all_edges_visited(self, edge_count):
        """
            Check if all edges in the edge status dictionary have been visited.

            Args:
                edge_count (Dict[Tuple[str, str], int]): edge status dictionary.

            Returns:
                bool: True if all edges have been visited, False otherwise.
        """
        for edge, count in edge_count.items():
            if count > 0:
                return False
        return True

    def _find_all_eulerian_paths(self, graph, node, path, edge_count):
        """
            Find all Eulerian paths in the given De Bruijn graph.

            Args:
                graph (Dict[str, List[str]]): De Bruijn graph.
                node (str): current node being visited.
                path (List[str]): current path being constructed.
                edge_count (Dict[Tuple[str, str], bool]): edge status dictionary.

            Yields:
                Generator[List[str], None, None]: a generator yielding all eulerian paths.
        """
        if self._is_all_edges_visited(edge_count):
            yield path[:]

        for neighbor in graph[node]:
            if edge_count[(node, neighbor)] <= 0:
                continue
            edge_count[(node, neighbor)] -= 1
            path.append(neighbor)
            yield from self._find_all_eulerian_paths(graph, neighbor, path, edge_count)
            path.pop()
            edge_count[(node, neighbor)] += 1

    def find_eulerian_path(self, graph):
        """
            Find an Eulerian path and all Eulerian paths in the given De Bruijn graph.

            Args:
                graph (Dict[str, List[str]]): The De Bruijn graph.

            Returns:
                Tuple[List[str], List[List[str]]
        """
        hierholzer_eulerian_path = []
        all_eulerian_paths = []

        is_eulerian, start_node = self._eulerian_check(graph)
        print(f"Eulerian Check: [{is_eulerian}]")

        if is_eulerian:
            # Hierholzer
            graph_copy = copy.deepcopy(graph)
            hierholzer_eulerian_path = self._hierholzer(graph_copy, start_node)

            # Backtrack
            edge_count = self._init_edge_count(graph)
            all_eulerian_paths = list(self._find_all_eulerian_paths(graph, start_node, [start_node], edge_count))
            unique_eulerian_paths = []
            for path in all_eulerian_paths:
                if path not in unique_eulerian_paths:
                    unique_eulerian_paths.append(path)
            all_eulerian_paths = unique_eulerian_paths

        return hierholzer_eulerian_path, all_eulerian_paths

    def assemble(self, eulerian_paths):
        """
            Assemble the given Eulerian paths into sequences.

            Args:
                eulerian_paths (List[List[str]]): list of eulerian paths

            Returns:
                List[str]: list of assembled sequences.
        """
        assembled_sequences = set()
        for path in eulerian_paths:
            sequence = path[0]
            for node in path[1:]:
                sequence += node[-1]
            assembled_sequences.add(sequence)
        return assembled_sequences

    def build_graph(self, fragments):
        """
            Build a De Bruijn graph from the given fragments.

            Args:
                fragments (List[str]): A list of fragments.

            Returns:
                Dict[str, List[str]]: De Bruijn graph.
        """
        graph = {}
        for fragment in fragments:
            prefix = fragment[:-1]
            suffix = fragment[1:]
            if prefix not in graph:
                graph[prefix] = []
            if suffix not in graph:
                graph[suffix] = []
            graph[prefix].append(suffix)
        return graph

    def save_as_svg(self, graph, filename="de_bruijn_graph.svg"):
        pydot_graph = pydot.Dot("de_bruijn_graph", graph_type="digraph", bgcolor="white")
        for u, adj in graph.items():
            node_u = pydot.Node(u, label=u, shape="circle")
            pydot_graph.add_node(node_u)

            for v in adj:
                node_v = pydot.Node(v, label=v, shape="circle")
                pydot_graph.add_edge(pydot.Edge(node_u, node_v, label=u[0]+v))

        pydot_graph.write_svg('de_brujin_graph.svg')
        return pydot_graph

def main():
    """
        Entry Point
    """

    sequence = "GACTTACGTACT"
    k = 3

    print(LINE)
    print(f"Sequence         : {sequence}")
    print(f"Fragment Size (k): {k}")

    fragment_generator = FragmentGenerator()
    fragments = fragment_generator.generate_fragments(sequence, k)

    print(LINE)
    print("Fragments: (Generated from Sequence)")
    print(fragments)

    de_bruijn_assembler = DeBruijnAssembler()
    de_bruijn_graph = de_bruijn_assembler.build_graph(fragments)

    print(LINE)
    print("De Bruijn Graph: (Built from Fragments)")
    print(f"Number of Nodes: {len(de_bruijn_graph)}")
    for node, neighbors in de_bruijn_graph.items():
        print(f"{node}: {neighbors}")
    print(LINE)

    hierholzer_eulerian_path, all_eulerian_paths = de_bruijn_assembler.find_eulerian_path(de_bruijn_graph)

    print("Eulerian Path: (Hierholzer Algorithm)")
    print(hierholzer_eulerian_path)
    print(LINE)

    print("All Eulerian Paths: (Backtracking Algorithm)")
    for eulerian_path in all_eulerian_paths:
        print(eulerian_path)
    print(LINE)

    assembled_sequences = de_bruijn_assembler.assemble(all_eulerian_paths)
    print("Assembled Sequences: (From All Eulerian Paths)")
    for assmbled_sequence in assembled_sequences:
        print(assmbled_sequence)
    print(LINE)

    de_bruijn_assembler.save_as_svg(de_bruijn_graph)

if __name__ == "__main__":
    main()
