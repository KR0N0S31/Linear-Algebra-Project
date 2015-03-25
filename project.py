import copy
import math
from pprint import pprint

class Matrix:
    def __init__(self, n, m):
        self.row = n
        self.col = m
        self.arr = [[0 for x in range(n)] for x in range(m)]

    def solve_system_for_triangular_matrix(self, b):
        solutions = Vector(self.row)
        if self.arr[0][self.col - 1] == 0: #upper-triangular
            for i in range(self.row):
                row_sum = 0
                for j in range(self.col):
                    if i == 0 and j == 0:
                        solutions.arr[i] = b.arr[i] / self.arr[i][j]
                    else:
                        if j < i:
                            row_sum += solutions.arr[j] * self.arr[i][j]
                        elif j == i:
                            rhs = b.arr[i] - row_sum
                            new_solution = rhs / self.arr[i][j]
                            solutions.arr[i] = new_solution
        elif self.arr[self.row - 1][0] == 0: #lower-triangular
            for i in range(self.row - 1, -1, -1):
                row_sum = 0
                for j in range(self.col - 1, -1, -1):
                    if i == self.row - 1 and j == self.col - 1:
                        solutions.arr[i] = b.arr[i] / self.arr[i][j]
                    else:
                        if j > i:
                            row_sum += solutions.arr[j] * self.arr[i][j]
                        elif j == i:
                            rhs = b.arr[i] - row_sum
                            new_solution = rhs / self.arr[i][j]
                            solutions.arr[i] = new_solution
        for i in range(self.row):
            solutions.arr[i] = round(solutions.arr[i], 14)
        return solutions

    def init_identity(self):
        for x in range(self.row):
            for y in range(self.col):
                if (x == y):
                    self.arr[x][y] = 1
                else:
                    self.arr[x][y] = 0

    def transpose(self):
        result = Matrix(self.row, self.col)
        for i in range(self.row):
            for j in range(self.col):
                result.arr[j][i] = self.arr[i][j]
        return result

    def multiply(self, other):
        result = Matrix(self.row, other.col)
        for i in range(len(self.arr)):
            for j in range(len(other.arr[0])):
                for k in range(len(other.arr)):
                    result.arr[i][j] += self.arr[i][k] * other.arr[k][j]
        return result

    def  multiply_by_vector(self, vector):
        result = Vector(self.row)
        for i in range(self.row):
            for j in range(self.col):
                result.arr[i] += self.arr[i][j] * vector.arr[j]
        return result

    def multiply_by_scalar(self, scalar):
        for i in range(self.row):
            for j in range(self.col):
                self.arr[i][j] *= scalar

    def divide_by_scalar(self, scalar):
        for i in range(self.row):
            for j in range(self.col):
                self.arr[i][j] /= scalar

    def subtract(self, other):
        result = Matrix(self.row, self.col)
        for i in range(len(self.arr)):
            for j in range(len(other.arr)):
                result.arr[i][j] = self.arr[i][j] - other.arr[i][j]
        return result

    def inf_norm(self):
        norm = float("-inf")
        for i in range(self.row):
            curr_sum = 0
            for j in range(self.col):
                curr_sum += self.arr[i][j]
            norm = curr_sum if (curr_sum > norm) else norm
        return norm

    def get_vector_at(self, k):
        v = Vector(self.row)
        for i in range(self.row):
            v.arr[i] = matrix.arr[i][k]
        return v

    def round(self):
        for i in range(self.row):
            for j in range(self.col):
                self.arr[i][j] = round(self.arr[i][j], 14)
                if self.arr[i][j] == -0:
                    self.arr[i][j] = 0

    def print(self):
        for i in range(self.row):
            for j in range(self.col):        
                if i == 0 and j == 0:
                    print("       [[%10.5f" % (self.arr[i][j]), sep="", end=" ")
                elif j == 0:
                    print("\t[%10.5f" % (self.arr[i][j]), sep="", end=" ")
                elif j == self.col - 1:
                    if i == self.row - 1:
                        print("%12.5f]]" % (self.arr[i][j]), sep="", end=" ")
                    else:
                        print("%12.5f]" % (self.arr[i][j]), sep="")
                else:
                    print("%12.5f" % (self.arr[i][j]), end="")
        print()

class Vector:
    def __init__(self, n):
        self.arr = [0] * n

    def norm(self):
        sum_of_squares = 0
        for k in self.arr:
            sum_of_squares = sum_of_squares + k * k
        return math.sqrt(sum_of_squares)

    def multiply_self_by_transpose(self):
        result = Matrix(len(self.arr), len(self.arr))
        for i in range(len(self.arr)):
            for j in range(len(self.arr)):
                result.arr[i][j] = self.arr[i] * self.arr[j]
        return result

def lu_fact(matrix):
    L = Matrix(matrix.row, matrix.col)
    L.init_identity()
    U = Matrix(matrix.row, matrix.col)
    U = copy.deepcopy(matrix)
    prev_row = 0
    curr_row = 1
    for col in range(U.col):
        for row in range(curr_row, U.row):
            if U.arr[row][col] != 0:
                scalar = U.arr[row][col] / U.arr[prev_row][col]
                L.arr[row][col] = scalar
                scalar *= -1
                for y in range(col, U.col):
                    first_entry = U.arr[prev_row][y]
                    victim_entry = U.arr[row][y]
                    new_value = first_entry * scalar + victim_entry
                    U.arr[row][y] = new_value
        prev_row = curr_row
        curr_row += 1
    print("----------------------LU Decomposition----------------------")
    print("A = ")
    matrix.print()
    print("L =")
    L.print()
    print("U =")
    U.print()
    print("L * U =")
    L.multiply(U).print()
    error_matrix = L.multiply(U).subtract(matrix)
    print("\nError = %f" % error_matrix.inf_norm())
    print("------------------------------------------------------------\n")
    return L, U

def qr_fact_househ(matrix):
    Q = Matrix(matrix.row, matrix.col)
    Q.init_identity()
    R = copy.deepcopy(matrix)
    for k in range(matrix.row - 1):
        a_i = R.get_vector_at(k)
        if (k > 0):
            for row_element in range(k):
                a_i.arr.pop(0)
        u = a_i
        u.arr[0] = a_i.arr[0] + a_i.norm()
        bottom = u.norm() * u.norm()
        top = u.multiply_self_by_transpose()
        top.multiply_by_scalar(2)
        top.divide_by_scalar(bottom)
        I = Matrix(matrix.row - k, matrix.col - k)
        I.init_identity()
        Hi = I.subtract(top)
        new_Hi = Matrix(matrix.row, matrix.col)
        for i in range(matrix.row):
            for j in range(matrix.col):
                if i > k - 1 and j > k - 1:
                    new_Hi.arr[i][j] = Hi.arr[i - k][j - k]
                else:
                    if i == j:
                        new_Hi.arr[i][j] = 1
                    else:
                        new_Hi.arr[i][j] = 0
        R = new_Hi.multiply(R)
        Q = Q.multiply(new_Hi)
        Q.round()
        R.round()
    print("-----------------------Householder QR-----------------------")
    print("A =")
    matrix.print()
    print("Q =")
    Q.print()
    print("R =")
    R.print()
    print("Q * R =")
    Q.multiply(R).print()
    error_matrix = Q.multiply(R).subtract(matrix)
    print("\nError = %f" % error_matrix.inf_norm())
    print("------------------------------------------------------------\n")
    return Q, R

def qr_fact_givens(matrix):
    Q = Matrix(matrix.row, matrix.col)
    Q.init_identity()
    R = copy.deepcopy(matrix)

    k = 1
    count = 0
    for j in range(matrix.col - 1):
        for i in range(k, matrix.row):
            top_entry = R.arr[k - 1][j]
            bottom_entry = R.arr[i][j]
            denom = math.sqrt(top_entry ** 2 + bottom_entry ** 2)
            c = top_entry / denom
            s = bottom_entry / denom
            Gi = Matrix(matrix.row, matrix.col)
            Gi.init_identity()
            Gi.arr[k - 1][j] = c
            Gi.arr[j][i] = s
            Gi.arr[i][j] = -s
            Gi.arr[i][i] = c
            R = Gi.multiply(R)
            R.round()
            Q = Q.multiply(Gi.transpose())
            Q.round()
        k += 1
    print("------------------------Givens QR---------------------------")
    print("A =")
    matrix.print()
    print("Q =")
    Q.print()
    print("R =")
    R.print()
    print("Q * R =")
    Q.multiply(R).print()
    error_matrix = Q.multiply(R).subtract(matrix)
    print("\nError = %f" % error_matrix.inf_norm())
    print("------------------------------------------------------------\n")
    return Q, R
    
def solve_lu_b(original, L, U, b):
    L_sol = L.solve_system_for_triangular_matrix(b)
    U_sol = U.solve_system_for_triangular_matrix(L_sol)
    print("------------------------------------------------------------")
    print("Using LU factorization to solve Ax = b for x:")
    print("A =")
    original.print()
    print("b =")
    print("\t", end="")
    pprint(b.arr)
    print("Result:")
    print("Ly = b:")
    print(" y = ", end="")
    pprint(L_sol.arr)
    print("     via forward substitution")
    print("\nUx = y:")
    print(" x = ", end="")
    pprint(U_sol.arr)
    print("     via backward substitution")
    print("------------------------------------------------------------")
    print()

def solve_qr_b(original, Q, R, b):
    Q_transpose_b = Q.transpose().multiply_by_vector(b)
    final_sol = R.solve_system_for_triangular_matrix(Q_transpose_b)
    print("------------------------------------------------------------")
    print("QR factorization to solve Ax = b for x:")
    print("A =")
    original.print()
    print("b =")
    print("\t", end="")
    pprint(b.arr)
    print("\nResult:\nRx = Q^t b")
    print("Rx = ", end="")
    pprint(Q_transpose_b.arr)
    print("x  = ", end="")
    pprint(final_sol.arr)
    print("     via backward substitution")
    print("------------------------------------------------------------")
    print()


matrix = Matrix(4, 4)
matrix.arr = [[1, 1, -1, 5], [1, -2, 3, 7], [2, 3, 1, 8], [5, 3, -1, 15]]
LU = lu_fact(matrix)
QR = qr_fact_givens(matrix)
b = Vector(matrix.row)
b.arr = [4, -6, 7, 77]
solve_lu_b(matrix, LU[0], LU[1], b)
solve_qr_b(matrix, QR[0], QR[1], b)


