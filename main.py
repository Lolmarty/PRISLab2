import copy
import sys
import string


class IntervalNumber:
    def __init__(self, low, high):
        self.low = low
        self.high = high

    def __str__(self):
        return "(" + str(self.low) + ";" + str(self.high) + ")"

    def Inverse(self):
        return IntervalNumber(1. / self.high, 1. / self.low)

    def DistanceToZero(self):
        return (((self.low + self.high) / 2.) ** 2 + (((self.high - self.low) / 2.) ** 2) / 3) ** 0.5


class TrapezoidalNumber:
    def __init__(self, low, mid_low, mid_high, high):
        self.low = low
        self.mid_low = mid_low
        self.mid_high = mid_high
        self.high = high

    def __str__(self):
        return "(" + str(self.low) + ";" + str(self.mid_low) + ";" + str(self.mid_high) + ";" + str(self.high) + ")"

    def AlphaLevel(self, alpha):
        return IntervalNumber(self.low + alpha * (self.mid_low - self.low),
                              self.high + alpha * (self.mid_high - self.high))

    def Inverse(self):
        return TrapezoidalNumber(
            1. / self.high,
            1. / self.mid_high,
            1. / self.mid_low,
            1. / self.low)


_1 = TrapezoidalNumber(1, 1, 1, 1)
one = TrapezoidalNumber(1. / 3, 0.5, 1.5, 2)
two = TrapezoidalNumber(1, 1.5, 2.5, 3)
three = TrapezoidalNumber(2, 2.5, 3.5, 4)
four = TrapezoidalNumber(3, 3.5, 4.5, 5)
five = TrapezoidalNumber(4, 4.5, 5.5, 6)
six = TrapezoidalNumber(5, 5.5, 6.5, 7)
seven = TrapezoidalNumber(6, 6.5, 7.5, 8)
eight = TrapezoidalNumber(7, 7.5, 8.5, 9)
nein = TrapezoidalNumber(8, 9, 9, 9)


class FuzzyPairwiseComparisonMatrix:
    def __init__(self, n, matrix):
        self.n = n
        self.matrix = copy.deepcopy(matrix)

    def __getitem__(self, i, j):
        return self.matrix[i][j]

    def __str__(self):
        return "\n".join(["[" + ",".join([str(element) for element in line]) + "]" for line in self.matrix])

    def AlphaLevel(self, alpha):
        return FuzzyPairwiseComparisonMatrix(self.n,
                                             [[self.matrix[i][j].AlphaLevel(alpha)  # if i != j else _1
                                               for j in range(0, self.n)]
                                              for i in range(0, self.n)])

    def GenerateFromRow(self, row_index):
        return FuzzyPairwiseComparisonMatrix(self.n,
                                             [[IntervalNumber(
                                                 self.matrix[i][row_index].low * self.matrix[row_index][j].low,
                                                 self.matrix[i][row_index].high * self.matrix[row_index][j].high)
                                               for j in range(0, self.n)]
                                              for i in range(0, self.n)])

    def Consistency(self):
        consistent = True
        for i in range(0, self.n):
            for j in range(0, self.n):
                if (i != j):
                    max_low_produce = max([self.matrix[i][k].low * self.matrix[k][j].low
                                           # if k != i and k != j else -sys.maxint - 1
                                           for k in range(0, self.n)])
                    min_high_produce = min([self.matrix[i][k].high * self.matrix[k][j].high
                                            # if k != i and k != j else sys.maxint
                                            for k in range(0, self.n)])
                    consistent = consistent and max_low_produce <= min_high_produce
        return consistent


class FuzzyWeights:
    def __init__(self, n, weights):
        """
        :type n: int
        :type weights: list[IntervalNumber]
        """
        if n != len(weights):
            raise ValueError("Length of the weights vector ({0}) doesn't match the number of dimensions ({1}).".format(
                len(weights), n))
        self.n = n
        self.weights = copy.deepcopy(weights)

    def GenerateSpectreElements(self):
        """Generates the h'th elements for spectres,
        transpose after finding them all"""
        return [weight.DistanceToZero() for weight in self.weights]


class FuzzyConsistencyCoefficientGenerator:
    def __init__(self, n, weights):
        """
        :type n: int
        :type weights: list[FuzzyWeights]
        """
        if n != len(weights):
            raise ValueError("Amount of the weights vectors ({0}) doesn't match the number of dimensions ({1}).".format(
                len(weights), n))
        if not all([n == weight.n for weight in weights]):
            indexes = ",".join([str(weights.index(weight)) for weight in weights if n != weight.n])
            raise ValueError("Length of the weights vectors #{0} doesn't match the number of dimensions ({1}).".format(
                indexes, n))
        self.n = n
        self.weights = copy.deepcopy(weights)

    def blarg(self):
        temp = [weight.GenerateSpectreElements() for weight in self.weights]
        spectres = [[temp[j][i] for j in range(0, self.n)] for i in range(0, self.n)]
        phis = []


criteria_fpcm = FuzzyPairwiseComparisonMatrix(3, [[_1, two.Inverse(), three.Inverse()],
                                                  [two, _1, three.Inverse()],
                                                  [three, three, _1]])

alternative_fpcm_by_crit_1 = FuzzyPairwiseComparisonMatrix(4, [[_1, one, one, three],
                                                               [one.Inverse(), _1, three, one],
                                                               [one.Inverse(), three.Inverse(), _1, two],
                                                               [three.Inverse(), one.Inverse(), two.Inverse(), _1]])

alternative_fpcm_by_crit_2 = FuzzyPairwiseComparisonMatrix(4, [[_1, three, three, three],
                                                               [three.Inverse(), _1, three.Inverse(), three.Inverse()],
                                                               [three.Inverse(), three, _1, one],
                                                               [three.Inverse(), three, one.Inverse(), _1]])

alternative_fpcm_by_crit_3 = FuzzyPairwiseComparisonMatrix(4, [[_1, one, three.Inverse(), two],
                                                               [one.Inverse(), _1, three.Inverse(), three],
                                                               [three, three, _1, five],
                                                               [two.Inverse(), three.Inverse(), five.Inverse(), _1]])

# for key in globals().keys():
#     if "fpcm" in key:
#         print key
#         print globals()[key]
#         print "alpha 0 "
#         print globals()[key].AlphaLevel(0)
#         print "consistent " + str(globals()[key].AlphaLevel(0).Consistency())
#         for i in range(0, globals()[key].n):
#             print "generated from row " + str(i)
#             print globals()[key].AlphaLevel(0).GenerateFromRow(i)
#         print "alpha 0.5 "
#         print globals()[key].AlphaLevel(0.5)
#         print "consistent " + str(globals()[key].AlphaLevel(0.5).Consistency())
#         for i in range(0, globals()[key].n):
#             print "generated from row " + str(i)
#             print globals()[key].AlphaLevel(0.5).GenerateFromRow(i)
#         print

dummy = FuzzyWeights(3, [IntervalNumber(0.1, 0.3), IntervalNumber(0.3, 0.4), IntervalNumber(0.2, 0.6)])
print dummy.GenerateSpectreElements()
