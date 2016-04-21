from __future__ import division
from numpy import *
from scipy.optimize import linprog

eps = 1e-5


class SimplexSolver:
    def __init__(self, obj):
        self.initial_vars_count = len(obj)
        self.objective = [1] + obj
        self.rows = []
        self.cons = []

    def add_constraint(self, expression, value):
        self.rows.append([0] + expression)
        self.cons.append(value)

    def _select_variable_for_pivot(self):
        low = 0
        idx = 0
        for i in range(1, len(self.objective) - 1):
            if self.objective[i] < low:
                low = self.objective[i]
                idx = i
        if idx == 0: return -1
        return idx

    def _select_row_for_pivot(self, col):
        right_hand_side = [self.rows[i][-1] for i in range(len(self.rows))]
        left_hand_side = [self.rows[i][col] for i in range(len(self.rows))]
        ratio = []
        for i in range(len(right_hand_side)):
            if abs(left_hand_side[i]) < eps:
                ratio.append(abs(max(right_hand_side)) / eps)
            else:
                ratio.append(right_hand_side[i] / left_hand_side[i])
        return argmin(ratio)

    def display(self):
        print '\n', matrix([self.objective] + self.rows)

    def _pivot(self, row, col):
        e = self.rows[row][col]
        self.rows[row] /= e
        for r in range(len(self.rows)):
            if r != row:
                self.rows[r] = self.rows[r] - self.rows[r][col] * self.rows[row]
        self.objective = self.objective - self.objective[col] * self.rows[row]

    def _check(self):
        if min(self.objective[1:-1]) >= 0: return 1
        return 0

    def _get_solutions(self):
        basic_vars_indexes = [index for index in range(1, len(self.objective) - 1)
                              if abs(self.objective[index]) < eps]  # these are the values that we are looking for

        basic_vars_rows = [[row[index] for index in basic_vars_indexes] for row in self.rows]
        transposed_basic_vars_rows = map(list, zip(*basic_vars_rows))
        basic_vars_value_indexes = [next(index for index in range(len(column)) if column[index] == 1)
                                    for column in transposed_basic_vars_rows]
        solutions = [0] * len(self.objective)
        for i in range(len(basic_vars_indexes)):
            solutions[basic_vars_indexes[i]] = self.rows[basic_vars_value_indexes[i]][-1]
        return solutions

    def prepare(self):
        # build full tableau
        for i in range(len(self.rows)):
            self.objective += [0]
            ident = [0 for r in range(len(self.rows))]
            ident[i] = 1
            self.rows[i] += ident + [self.cons[i]]
            self.rows[i] = array(self.rows[i], dtype=float)
        self.objective = array(self.objective + [0], dtype=float)

    def solve(self):
        # solve
        # self.display()
        while not self._check():
            c = self._select_variable_for_pivot()
            r = self._select_row_for_pivot(c)
            self._pivot(r, c)
            # print '\npivot column: %s\npivot row: %s' % (c + 1, r + 2)
            # self.display()
        solutions = self._get_solutions()
        return solutions[1:self.initial_vars_count + 1]


if __name__ == '__main__':
    """
    max f = 12x + 40y
    st
    x + y <= 16
    x + 3y <= 36
    x <= 10
    x,y >= 0
    """

    t = SimplexSolver([-12, -40])
    t.add_constraint([1, 1], 16)
    t.add_constraint([1, 3], 36)
    t.add_constraint([1, 0], 10)
    print t.solve()


def interval_weights(A):
    n = len(A)
    c = [1] * n * (n - 1)
    c += [0] * n
    B = [0] * n * (n - 1)
    b = [0] * n * (n - 1)
    for i in range(0, n * (n - 1)):
        B[i] = [0] * n * n
    l = 0
    for i in range(0, n * (n - 1)):
        B[i][i] = -1
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            B[l][n * (n - 1) + i] = -1
            B[l][n * (n - 1) + j] = 1
            B[round(l + n * (n - 1) / 2)][n * (n - 1) + i] = 1
            B[round(l + n * (n - 1) / 2)][n * (n - 1) + j] = -1
            b[l] = - math.log(A[i][j][0])
            b[round(l + n * (n - 1) / 2)] = math.log(A[i][j][1])
            l += 1
    l = [0] * n * (n - 1)
    l += [None] * n
    u = [None] * n * (n - 1)
    u += [0] * n
    delta = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n * n)]))
    l = 0
    print('Before:')
    disp_inter_matr(A)
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            A[i][j][0] *= math.e ** (- delta.x[l])
            A[i][j][1] *= math.e ** delta.x[round(l + n * (n - 1) / 2)]
            l += 1
    print('After:')
    disp_inter_matr(A)
    B = [0] * n * (n - 1)
    for i in range(0, n * (n - 1)):
        B[i] = [0] * n
    l = 0
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            B[l][i] = -1
            B[l][j] = 1
            B[round(l + n * (n - 1) / 2)][i] = 1
            B[round(l + n * (n - 1) / 2)][j] = -1
            b[l] = - math.log(A[i][j][0])
            b[round(l + n * (n - 1) / 2)] = math.log(A[i][j][1])
            l += 1
    l = [None] * n
    u = [0] * n
    w = [0] * n
    for i in range(0, n):
        w[i] = [0] * 2
        c = [0] * n
        c[i] = 1
        t = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n)]))
        tl = [0] * n
        for j in range(0, n):
            tl[j] = math.e ** t.x[j]
        tl = [x / sum(tl) for x in tl]
        w[i][0] = tl[i]
        c[i] = -1
        t = linprog(c=c, A_ub=B, b_ub=b, bounds=tuple([(l[i], u[i]) for i in range(0, n)]))
        tl = [0] * n
        for j in range(0, n):
            tl[j] = math.e ** t.x[j]
        tl = [x / sum(tl) for x in tl]
        w[i][1] = tl[i]
    return w


def disp_inter_matr(A):
    for i in A:
        for j in i:
            print('[%0.2f; %0.2f] ' % (j[0], j[1]))
        print()
