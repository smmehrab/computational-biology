# ************************************************
# username  :   smmehrab
# fullname  :   s.m.mehrabul islam
# email     :   smmehrabul-2017614964@cs.du.ac.bd
# institute :   university of dhaka, bangladesh
# reg       :   2017614964
# ************************************************

from typing import Tuple

DEBUG = True
LINE = "-----------------------------------------"

class Waterman:
    """
        Smith-Waterman Algorithm
        for Pairwise Local Alignment
        (Dynamic Programming)

        Methods:
            - __init__ (constructor)
            - set_sequences
            - set_scoring_parameters
            - find_optimal_alignment
    """

    def __init__(self, sequence1: str, sequence2: str, match: int, mismatch: int, gap: int):
        """
            Constructor
        """
        self.set_sequences(sequence1, sequence2)
        self.set_scoring_parameters(match, mismatch, gap)

    def set_sequences(self, sequence1: str, sequence2: str):
        """
            Set the sequences
        """
        self.sequence1 = sequence1
        self.sequence2 = sequence2
        self.sequence1_length = len(self.sequence1)
        self.sequence2_length = len(self.sequence2)
        # matrix init
        self.matrix = []
        for i in range(self.sequence1_length+1):
            self.matrix.append([])
            for _ in range(self.sequence2_length+1):
                self.matrix[i].append(0)

    def set_scoring_parameters(self, match: int, mismatch: int, gap: int):
        """
            Set the scoring parameters
        """
        if not match or not mismatch or not gap:
            raise Exception("Invalid Scoring Parameters")
        self.scores = {
            "match": match,
            "mismatch": mismatch,
            "gap": gap
        }

    def _set_gap_row_column(self):
        """
            Set the gap row and column
        """
        for i in range(self.sequence1_length+1):
            self.matrix[i][0] = 0
        for j in range(self.sequence2_length+1):
            self.matrix[0][j] = 0

    def _populate_matrix(self):
        """
            Populate the matrix (DP)
        """
        # gap row and gap column
        self._set_gap_row_column()
        # rest of the matrix
        for i in range(1, self.sequence1_length+1):
            for j in range(1, self.sequence2_length+1):
                # increment both sequences
                current_score = (self.scores["match"] if self.sequence1[i-1] == self.sequence2[j-1] else self.scores["mismatch"])
                diagonal_score = self.matrix[i-1][j-1] + current_score
                # increment only sequence1 and gap in sequence2
                vertical_score = self.matrix[i-1][j] + self.scores["gap"]
                # increment only sequence2 and gap in sequence1
                horizontal_score = self.matrix[i][j-1] + self.scores["gap"]
                self.matrix[i][j] = max(0, diagonal_score, horizontal_score, vertical_score)
        if DEBUG:
            self._print_matrix()

    def _find_alignment_score(self, aligned_sequence1: str, aligned_sequence2: str) -> int:
        """
            Find the alignment score
            from two aligned sequences
        """
        alignment_score = 0
        alignment_length = len(aligned_sequence1)
        for i in range(alignment_length):
            if aligned_sequence1[i] == aligned_sequence2[i]:
                alignment_score += self.scores["match"]
            elif aligned_sequence1[i] == '-' or aligned_sequence2[i] == '-':
                alignment_score += self.scores["gap"]
            else:
                alignment_score += self.scores["mismatch"]
        return alignment_score

    def _traceback(self) -> Tuple[str, str, int]:
        """ 
            Traceback to find the optimal alignment
            from the matrix
        """
        aligned_sequence1 = ""
        aligned_sequence2 = ""

        if DEBUG:
            print(LINE)
            print("Traceback")
            print(LINE)

        # Find the max score
        max_score = 0
        max_i, max_j = 0, 0
        for i in range(self.sequence1_length + 1):
            for j in range(self.sequence2_length + 1):
                if self.matrix[i][j] > max_score:
                    max_score = self.matrix[i][j]
                    max_i, max_j = i, j

        i, j = max_i, max_j
        while i > 0 and j > 0 and self.matrix[i][j] > 0:
            if DEBUG:
                print(f"({self.matrix[i][j]}) --> ", end="")
            # diagonal
            diagonal_score = self.matrix[i - 1][j - 1] + (self.scores["match"] if self.sequence1[i - 1] == self.sequence2[j - 1] else self.scores["mismatch"])
            if self.matrix[i][j] == diagonal_score:
                aligned_sequence1 = self.sequence1[i - 1] + aligned_sequence1
                aligned_sequence2 = self.sequence2[j - 1] + aligned_sequence2
                i -= 1
                j -= 1
            # vertical
            elif self.matrix[i - 1][j] + self.scores["gap"] == self.matrix[i][j]:
                aligned_sequence1 = self.sequence1[i - 1] + aligned_sequence1
                aligned_sequence2 = '-' + aligned_sequence2
                i -= 1
            # horizontal
            else:
                aligned_sequence1 = '-' + aligned_sequence1
                aligned_sequence2 = self.sequence2[j - 1] + aligned_sequence2
                j -= 1
        if DEBUG:
            print(f"({self.matrix[i][j]})")

        alignment_score = max_score
        return aligned_sequence1, aligned_sequence2, alignment_score

    def _print_matrix(self):
        """
            Print the matrix (Formatted)
        """
        # title
        print(LINE)
        print("Matrix (vertical: sequence1, horizontal: sequence2))")
        print(LINE)
        # header row
        box = "□"
        print(f"{box:>5} {'-':>5}", end=" ")
        for c in self.sequence2:
            print(f"{c:>5}", end=" ")
        print()
        # normal rows
        for i in range(len(self.sequence1) + 1):
            # header column
            if i == 0:
                c = '-'
            # normal column
            else:
                c = self.sequence1[i - 1]
            print(f"{c:>5}", end=" ")
            for j in range(len(self.sequence2) + 1):
                # matrix value
                print(f"{self.matrix[i][j]:>5}", end=" ")
            print()

    def find_optimal_alignment(self) -> Tuple[str, str, int]:
        """
            Find Optimal Alignment
            Using Dynamic Programming
        """
        if not self.sequence1 or not self.sequence2:
            raise Exception("Sequences not set")
        if not self.scores:
            raise Exception("Scoring parameters not set")

        self._populate_matrix()
        return self._traceback()

def main():
    """
        Entry Point for the Program
    """
    
    # Input
    S1 = "CGTGAATTCAT"
    S2 = "GACTTAC"
    match = 5
    mismatch = -3
    gap = -4

    # Solve
    waterman = Waterman(S1, S2, match, mismatch, gap)
    aligned_sequence1, aligned_sequence2, score = waterman.find_optimal_alignment()

    # Output
    print(LINE)
    print(f"Alignment (Score: {score})")
    print(LINE)
    print(f"Sequence 1: {S1}")
    print(f"Sequence 2: {S2}")
    print(f"Locally Aligned Sequence 1: {aligned_sequence1}")
    print(f"Locally Aligned Sequence 2: {aligned_sequence2}")
    print(LINE)

if __name__ == "__main__":
    main()
