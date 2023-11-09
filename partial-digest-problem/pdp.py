# ************************************************
# username  :   smmehrab
# fullname  :   s.m.mehrabul islam
# email     :   smmehrabul-2017614964@cs.du.ac.bd
# institute :   university of dhaka, bangladesh
# reg       :   2017614964
# ************************************************

"""
    Write and algorithm to solve partial digest problem. Find X such that ΔX = L
    For example L = {2, 2, 3, 3, 4, 5, 6, 7, 8, 10}. Solve for L (i.e. find
    X such that ΔX = L).
"""

class PartialDigestProblem:

    def __init__(self):
        self.width = 0

    def _delta(self, target, numbers):
        """
            Finds the difference between target and each number in numbers
        """
        return [abs(target - number) for number in numbers]

    def _place(self, L, X):
        """
            Places the next number from L to X
            and recursively calls itself until
            a valid solution is found
        """
        # solution found
        if not L:
            print(X)
            return

        # 2 candidates = 2 branches
        candidates = [max(L), self.width - max(L)]
        for y in candidates:
            delta_yX = self._delta(y, X)
            if set(delta_yX).issubset(set(L)):
                for member in delta_yX:
                    L.remove(member)
                X.append(y)
                self._place(L, X)
                X.remove(y)
                L.extend(delta_yX)

    def solve(self, L):
        """
            Entry point to the algorithm
        """
        print(f"L: {L}")
        print("X: (Valid Solutions)")
        self.width = max(L)
        L.remove(self.width)
        X = [0, self.width]
        self._place(L, X)

if __name__ == "__main__":
    L = [2, 2, 3, 3, 4, 5, 6, 7, 8, 10]
    pdp = PartialDigestProblem()
    pdp.solve(L)
