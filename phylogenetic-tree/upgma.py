# ************************************************
# username  :   smmehrab
# fullname  :   s.m.mehrabul islam
# email     :   smmehrabul-2017614964@cs.du.ac.bd
# institute :   university of dhaka, bangladesh
# reg       :   2017614964
# ************************************************

# For Graph Visualization
# - sudo apt install graphviz
# - pip install pydot

import pydot
from typing import List, Tuple

DEBUG = True
LINE = "-----------------------------------------"

class UPGMA:
    """
        UPGMA Algorithm
        Methods:
            - __init__ (constructor)
            - set_sequences
            - build_phylogenetic_tree
            - save_tree
    """

    def __init__(self, sequences: List[str]):
        """
            Constructor
        """
        self.set_sequences(sequences)

    def set_sequences(self, sequences: List[str]):
        """
            Initialize the sequences
        """
        self.sequences = sequences
        self.number_of_sequence = len(sequences)
        self.sequence_length = len(sequences[0])

        self.distances = {}
        self.components = []
        self.tree = []

        self._populate_distance_matrix()

    def build_phylogenetic_tree(self) -> List[Tuple[str, str, float]]:
        """
            Build the phylogenetic tree
        """
        loop_count = self.number_of_sequence
        while loop_count > 1:
            a, b, distance = self._get_min_distance()
            self.tree.append((a, b, distance))
            self._update_distance_matrix(a, b)
            loop_count -= 1
        print(LINE)
        return self.tree

    def _populate_distance_matrix(self):
        """
            Populate the distance matrix
        """
        for i in range(self.number_of_sequence):
            si = str(i+1)
            self.components.append(si)
            for j in range(i+1, self.number_of_sequence):
                sj = str(j+1)
                self.distances[(si, sj)] = self._get_distance(i, j)

    def _get_distance(self, a: int, b: int) -> int:
        """
            Get the distance between two sequences
            a: sequence index
            b: sequence index
        """
        distance = 0
        for k in range(self.sequence_length):
            if self.sequences[a][k] != self.sequences[b][k]:
                distance += 1
        return distance

    def _get_min_distance(self) -> Tuple[str, str, float]:
        """
            Get the minimum distance from the distance matrix
            a: sequence index
            b: sequence index
        """
        a, b, distance = None, None, float("inf")
        for pair, d in self.distances.items():
            if d < distance:
                distance = d
                a, b = pair
        return a, b, distance

    def _update_distance_matrix(self, a: str, b: str):
        """
            Update distance matrix after merging (connecting) two sequences
            a: sequence index
            b: sequence index
        """
        if DEBUG:
            print(f"Merging {a} and {b}:")

        merged = ''.join(sorted(a+b))
        updated_distances = {}
        for key, value in self.distances.items():
            if a not in key and b not in key:
                updated_distances[key] = value
        # Remove Old Components
        self.components.remove(a)
        self.components.remove(b)
        for component in self.components:
            # Keys
            key_a = (a, component) if (a, component) in self.distances else (component, a)
            key_b = (b, component) if (b, component) in self.distances else (component, b)
            # Distances (Avg)
            distance_a = self.distances[key_a] if key_a in self.distances else 0
            distance_b = self.distances[key_b] if key_b in self.distances else 0
            avg_distance = (distance_a + distance_b) / 2.00
            # Update
            updated_distances[(merged, component)] = avg_distance

            if DEBUG:
                print(f"    - Between {key_a}: {distance_a} and {key_b}: {distance_b}")
                print(f"    - Avg: {updated_distances[(merged, component)]}")

        # Add Merged Component
        self.components.append(merged)
        # Update Distance Matrix
        self.distances = updated_distances

        if DEBUG:
            print(f"    - Components: {self.components}")
            print(f"    - Distances: {self.distances}")

    def _get_parent_node(self, graph, parent):
        """
            Pydot
            Get parent node from the graph
        """
        for node in graph.get_nodes():
            if node.get_name() == parent:
                return node
        return pydot.Node(parent, label=parent, shape='point')

    def save_tree(self, filename: str = "phylogenetic_tree.png"):
        """
            Pydot
            Save a phylogenetic tree image as PNG
        """
        # Pydot Graph Config
        graph = pydot.Dot("phylogenetic_tree", graph_type="graph", bgcolor="white", rankdir="TB")
        graph.set_node_defaults(fillcolor="lightblue", style="filled", shape="circle", fontname="Courier", fontsize="10")
        graph.set_edge_defaults(fontname="Courier", fontsize="10")
        for a, b, _ in self.tree:
            parent = ''.join(sorted(a+b))
            # Nodes
            node_a = pydot.Node(a, label=a)
            node_b = pydot.Node(b, label=b)
            node_parent = self._get_parent_node(graph, parent)
            # Edges
            parent_to_a = pydot.Edge(node_parent, node_a)
            parent_to_b = pydot.Edge(node_parent, node_b)
            # Update Tree
            graph.add_node(node_a)
            graph.add_node(node_b)
            graph.add_node(node_parent)
            graph.add_edge(parent_to_a)
            graph.add_edge(parent_to_b)
        # Save
        graph.write_png(filename)

def main():
    """
        Entry Point for the Program
    """

    S1 = "ACGCGTTGGGCGATGGCAAC"
    S2 = "ACGCGTTGGGCGACGGTAAT"
    S3 = "ACGCATTGAATGATGATAAT"
    S4 = "ACGCATTGAATGATGATAAT"
    S5 = "ACACATTGAGTGATAATAAT"
    sequences = [S1, S2, S3, S4, S5]

    upgma = UPGMA(sequences)
    tree = upgma.build_phylogenetic_tree()

    print("Phylogenetic Tree:")
    print(LINE)
    print("{:<15} {:<15} {:}".format("node1", "node2", "distance"))
    print(LINE)
    for a, b, distance in tree:
        print("{:<15} {:<15} {:}".format(a, b, distance))
    print(LINE)

    upgma.save_tree()

if __name__ == "__main__":
    main()
