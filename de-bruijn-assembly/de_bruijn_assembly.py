# ************************************************
# username  :   smmehrab
# fullname  :   s.m.mehrabul islam
# email     :   smmehrabul-2017614964@cs.du.ac.bd
# institute :   university of dhaka, bangladesh
# reg       :   2017614964
# ************************************************

import copy
import pydot

class FragmentGenerator:
    def __init__(self, sequence=None, k=None):
        self._sequence = sequence
        self._k = k

    def set(self, sequence, k):
        """
            Set the sequence and fragment size (k)
        """
        self._sequence = sequence
        self._k = k

    def generate(self) :
        """
            Generate fragments according to the sequence and fragment size (k)
        """
        if self._sequence is None or self._k is None:
            raise Exception("Sequence or Fragment Size (k) is not set.")

        fragments = []
        sequence_length = len(self._sequence)
        for i in range(0, sequence_length - self._k + 1):
            fragment = self._sequence[i: i + self._k]
            fragments.append(fragment)
        return fragments

class DeBruijnAssembler:
    def __init__(self):
        pass

    def _find_degrees(self, graph):
        """
            Find indegree and outdegree of each node in the graph
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
            Check if the given De Bruijn Graph has an Eulerian path.
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
            Hierholzer: Find an Eulerian path in the given De Bruijn graph
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

    def _init_edge_status(self, graph):
        """
            Initialize the edge status for the given De Bruijn graph
            alongside getting the total number of edges
        """
        visited = {}
        total_edges = 0
        for node, neighbors in graph.items():
            visited[node] = []
            for _ in neighbors:
                total_edges += 1
                visited[node].append(False)
        return visited, total_edges

    def _find_all_eulerian_paths(self, graph, node, path, visited, edges_included, total_edges):
        """
            Backtrack: Find all Eulerian paths in the given De Bruijn graph
        """

        if edges_included == total_edges:
            yield path[:]

        for index, neighbor in enumerate(graph[node]):
            if not visited[node][index]:
                visited[node][index] = True
                path.append(neighbor)
                yield from self._find_all_eulerian_paths(graph, neighbor, path, visited, edges_included+1, total_edges)
                path.pop()
                visited[node][index] = False

    def find_eulerian_path(self, graph):
        """
            Find eulerian paths for the given De Bruijn graph
            (1) Eulerian Check
            (2) Hierholzer
            (3) Backtrack
        """

        # Eulerian Check
        is_eulerian, start_node = self._eulerian_check(graph)
        if not is_eulerian:
            raise Exception("Given De Bruijn Graph is not Eulerian.")

        # Hierholzer
        graph_copy = copy.deepcopy(graph)
        hierholzer_eulerian_path = self._hierholzer(graph_copy, start_node)

        # Backtrack
        visited, total_edges = self._init_edge_status(graph)
        all_eulerian_paths = list(self._find_all_eulerian_paths(graph, start_node, [start_node], visited, 0, total_edges))

        # Remove Duplicates
        unique_eulerian_paths = []
        for path in all_eulerian_paths:
            if path not in unique_eulerian_paths:
                unique_eulerian_paths.append(path)
        all_eulerian_paths = unique_eulerian_paths

        return hierholzer_eulerian_path, all_eulerian_paths

    def assemble(self, eulerian_paths):
        """
            Assemble the given Eulerian paths into sequences
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
            Build a De Bruijn graph from the given fragments
            and store it as a python dictionary
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

    def save_graph(self, graph, filename="de_bruijn_graph.png"):
        """
            Save a De Bruijn graph as an .png image.
        """
        # Graph Config
        pydot_graph = pydot.Dot("de_bruijn_graph", graph_type="digraph", bgcolor="white", rankdir="LR")
        pydot_graph.set_node_defaults(fillcolor="lime", style="filled", shape="circle", fontname="Courier", fontsize="10")
        pydot_graph.set_edge_defaults(fontname="Courier", fontsize="10")

        for node, neighbors in graph.items():
            # Node
            pydot_node = pydot.Node(node, label=node)
            pydot_graph.add_node(pydot_node)
            for neighbor in neighbors:
                # Neighbor
                pydot_neighbor = pydot.Node(neighbor, label=neighbor)
                edge_label = node[:-1]+neighbor
                # Edge
                pydot_edge = pydot.Edge(pydot_node, pydot_neighbor, label=edge_label)
                pydot_graph.add_edge(pydot_edge)
        # Save
        pydot_graph.write_png(filename)
        return pydot_graph

def output(sequence, k, fragments, graph, hierholzer_eulerian_path, all_eulerian_paths, assembled_sequences):
    LINE = "-------------------------------------------------"

    print(LINE)
    print(f"Sequence         : {sequence}")
    print(f"Fragment Size (k): {k}")

    print(LINE)
    print("Fragments: (Generated from Sequence)")
    print(fragments)

    print(LINE)
    print("De Bruijn Graph: (Built from Fragments)")
    print(f"Number of Nodes: {len(graph)}")
    for node, neighbors in graph.items():
        print(f"{node}: {neighbors}")
    print(LINE)

    print("Eulerian Path: (Hierholzer Algorithm)")
    print(hierholzer_eulerian_path)
    print(LINE)

    print("All Eulerian Paths: (Backtracking Algorithm)")
    for eulerian_path in all_eulerian_paths:
        print(eulerian_path)
    print(LINE)

    print("Assembled Sequences: (From All Eulerian Paths)")
    for assmbled_sequence in assembled_sequences:
        print(assmbled_sequence)
    print(LINE)

def main():
    """
        Entry Point for the Program
    """

    # Input
    sequence = "GACTTACGTACT"
    k = 3

    # Generate Fragments
    fragment_generator = FragmentGenerator(sequence, k)
    fragments = fragment_generator.generate()

    # Build Graph
    de_bruijn_assembler = DeBruijnAssembler()
    de_bruijn_graph = de_bruijn_assembler.build_graph(fragments)

    # Save Graph
    de_bruijn_assembler.save_graph(de_bruijn_graph)

    # Find Eulerian Paths
    hierholzer_eulerian_path, all_eulerian_paths = de_bruijn_assembler.find_eulerian_path(de_bruijn_graph)

    # Assemble into Sequences
    assembled_sequences = de_bruijn_assembler.assemble(all_eulerian_paths)

    # Print Output
    output(sequence, k, fragments, de_bruijn_graph, hierholzer_eulerian_path, all_eulerian_paths, assembled_sequences)

if __name__ == "__main__":
    main()
