class IntervalNumber:
    def __init__(self, low, high):
        self.low = low
        self.high = high

    def Inverse(self):
        return IntervalNumber(1. / self.high, 1. / self.low)


class TrapezoidalNumber:
    def __init__(self, low, mid_low, mid_high, high):
        self.low = low
        self.mid_low = mid_low
        self.mid_high = mid_high
        self.high = high

    def AlphaLevel(self, alpha):
        return IntervalNumber(self.low + alpha * (self.mid_low - self.low),
                              self.high + alpha * (self.mid_high - self.high))

    def Inverse(self):
        raise TrapezoidalNumber(
            1. / self.high,
            1. / self.mid_high,
            1. / self.mid_low,
            1. / self.low)


class FuzzyPairwiseComparisonMatrix:
    def __init__(self, n):
        self.n = n
        self.matrix = [[None for _ in range(0, n)] for __ in range(0, n)]

    def __getitem__(self, i, j):
        if self.matrix[i][j] is not None:
            return self.matrix[i][j]
        else:
            return self.matrix[j][i].Inverse()



one = TrapezoidalNumber(1./3, 0.5, 1.5, 2)
two = TrapezoidalNumber(1, 1.5, 2.5, 3)
three = TrapezoidalNumber(2, 2.5, 3.5, 4)
four = TrapezoidalNumber(3, 3.5, 4.5, 5)
five = TrapezoidalNumber(4, 4.5, 5.5, 6)
six = TrapezoidalNumber(5, 5.5, 6.5, 7)
seven = TrapezoidalNumber(6, 6.5, 7.5, 8)
eight = TrapezoidalNumber(7, 7.5, 8.5, 9)
nein = TrapezoidalNumber(8, 9, 9, 9)