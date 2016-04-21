import copy
import math
import simplex
import string

N = 11  # scale coefficient


class IntervalNumber:
    def __init__(self, low, high):
        self.low = low
        self.high = high

    def __str__(self):
        return "(" + str(self.low) + ";" + str(self.high) + ")"

    def __contains__(self, item):
        """
        :type item: float
        """
        return self.low <= item < self.high

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

    def GenerateMinimalExpandedMatrix(self):
        """
        delta_1 in the front, delta_2 in the middle, weights in the end
        :return:
        """
        solver = simplex.SimplexSolver([1] * (self.n * (self.n - 1)) + [0] * self.n)
        for i in range(self.n - 1):
            for j in range(i + 1, self.n):
                expression1 = [0] * self.n ** 2
                expression1[(self.n * 2 - 1 - i) * i / 2 + j] = 1
                expression1[self.n * (self.n - 1) + i] = -1
                expression1[self.n * (self.n - 1) + j] = 1
                value1 = math.log(self.matrix[i][j].low)
                solver.add_constraint(expression1, value1)

                expression2 = [0] * self.n ** 2
                expression2[self.n * (self.n - 1) / 2 + (self.n * 2 - 1 - i) * i / 2 + j] = 1
                expression2[self.n * (self.n - 1) + i] = 1
                expression2[self.n * (self.n - 1) + j] = 1
                value2 = -math.log(self.matrix[i][j].high)
                solver.add_constraint(expression2, value2)
        solver.prepare()
        solver.display()
        deltas_and_weights = solver.solve()
        print(deltas_and_weights)


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

    def GeneratePreSpectreElements(self):
        """Generates the h'th elements for spectres,
        transpose after finding them all"""
        return [weight.DistanceToZero() for weight in self.weights]


class Spectre:
    def __init__(self, n_scale, n_experts, pre_spectre):
        """
        :type n_scale: int
        :type n_experts: int
        :type pre_spectre: list[float]
        """
        if n_experts != len(pre_spectre):
            raise ValueError("Length of the pre spectre vector ({0}) doesn't match the number of experts ({1}).".format(
                len(pre_spectre), n_experts))
        self.n_experts = n_experts
        n_scale -= 1
        scale = [IntervalNumber(0., 0.5 / n_scale)]
        scale.extend([IntervalNumber(float(i) / n_scale - 0.5 / n_scale,
                                     float(i) / n_scale + 0.5 / n_scale)
                      for i in range(1, n_scale)])
        scale.append(IntervalNumber(1 - 0.5 / n_scale, 1))
        n_scale += 1
        self.n_scale = n_scale
        self.spectre = [0 for _ in range(0, n_scale)]
        for item in pre_spectre:
            item_index = [scale.index(interval) for interval in scale if item in interval][0]
            self.spectre[item_index] += 1

    def Phi(self):
        return -sum([math.log(float(item) / self.n_experts) * float(item) / self.n_experts
                     for item in self.spectre if item != 0])

    def Average(self):
        return float(sum([self.spectre[i] * i for i in range(0, self.n_scale)])) / self.n_experts

    def Psi(self):
        average = self.Average()
        return float(sum([self.spectre[i] * abs(i - average) for i in range(0, self.n_scale)])) / self.n_experts

    def HetaZero(self):
        G = self.n_experts / (self.n_scale * math.log(self.n_scale) * math.log(self.n_experts))
        return math.log(self.n_scale) + G * sum([abs(k - (self.n_scale + 1.) / 2.) for k in range(0, self.n_scale)])

    def Heta(self):
        return self.Phi() + self.Psi()

    def ConsistencyCoefficient(self):
        return 1 - self.Heta() / self.HetaZero()


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
        temp = [weight.GeneratePreSpectreElements() for weight in self.weights]
        pre_spectres = [[temp[j][i] for j in range(0, self.n)] for i in range(0, self.n)]
        spectres = [Spectre(N, self.n, pre_spectre) for pre_spectre in pre_spectres]


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
for matrix in [criteria_fpcm,
               alternative_fpcm_by_crit_1,
               alternative_fpcm_by_crit_2,
               alternative_fpcm_by_crit_3]:
    for alpha in [0., 0.5]:
        # print key
        # print globals()[key]
        # print "alpha " +str(alpha)
        # print globals()[key].AlphaLevel(alpha)
        # print "consistent " + str(globals()[key].AlphaLevel(alpha).Consistency())
        if not matrix.AlphaLevel(alpha).Consistency():
            matrix.AlphaLevel(alpha).GenerateMinimalExpandedMatrix()
            # for i in range(0, globals()[key].n):
            #     print "generated from row " + str(i)
            #     print globals()[key].AlphaLevel(alpha).GenerateFromRow(i)
    print

# dummy = FuzzyWeights(3, [IntervalNumber(0.1, 0.3), IntervalNumber(0.3, 0.4), IntervalNumber(0.2, 0.6)])
# print dummy.GeneratePreSpectreElements()

# spectre = Spectre(11, 5, [0.1, 0.499965, 0.7897, 0.8, 0.0000001])
# print spectre.spectre
# print spectre.Average()
# print spectre.Psi()
# print spectre.Phi()
# print spectre.HetaZero()
# print spectre.Heta()
# print spectre.ConsistencyCoefficient()
