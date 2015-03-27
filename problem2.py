import random
from pprint import pprint

class Matrix:
    def __init__(self, n, m):
        self.row = n
        self.col = m
        self.mat = [[0 for y in range(self.col)] for x in range(self.row)]
        self.b   = [0] * self.row

    def augment_with(self, b):
        self.b = b

class Vector:
    def __init__(self, n):
        self.row = n
        self.col = 1
        self.arr = [0] * n

class ConvolutionalMatrix:
    def __init__(self, n):
        self.n = n
        self.row = self.n + 3
        self.col = self.n

    def gen_random_x_stream(self):
        binary_stream = Vector(self.n)
        random.seed()
        for i in range(self.n):
            binary_stream.arr[i] = random.randint(0, 1)
        return binary_stream

    def gen_y_stream(self, x_stream):
        A0 = Matrix(self.row, self.col)
        A1 = Matrix(self.row, self.col)
        jth_pos = 0
        # Generate A0 and A1 simultaneously
        for x in range(self.row):
            for y in range(self.col):
                # Formula for generating A0
                if y == jth_pos or y == jth_pos - 2 or y == jth_pos - 3:
                    A0.mat[x][y] = 1
                else:
                    A0.mat[x][y] = 0
                # Formula for generating A1
                if y == jth_pos or y == jth_pos - 1 or y == jth_pos - 3:
                    A1.mat[x][y] = 1
                else:
                    A1.mat[x][y] = 0
            jth_pos += 1
        y0 = self.multiply_binary_matrix_by_vector(A0, x_stream)
        y1 = self.multiply_binary_matrix_by_vector(A1, x_stream)
        y  = self.merge_y0_with_y1(y0, y1)
        return y

    def multiply_binary_matrix_by_vector(self, matrix, vector):
        result = Vector(matrix.row)
        for i in range(matrix.row):
            for j in range(vector.row):
                result.arr[i] += matrix.mat[i][j] * vector.arr[j]
                result.arr[i] %= 2
        return result

    def merge_y0_with_y1(self, y0, y1):
        y = Vector(y0.row)
        for x in range(y0.row):
            concat = str(y0.arr[x]) + str(y1.arr[x])
            y.arr[x] = concat
        return y

def jacobi(matrix, initial_guess, tol):
    A = matrix.mat
    b = matrix.b
    iterations = Vector(matrix.row)
    xi = 0
    for k in range(20):
        for i in range(matrix.row):
            for j in range(matrix.col):
                if j == xi:
                    iterations.arr[i] += b[i] / A[i][xi]
                    iterations.arr[i] = round(iterations.arr[i], 14)
                    print(b[i], "/", A[i][xi], "=", iterations.arr[i], "=", b[i] / A[i][xi])
                else:
                    iterations.arr[i] += (-A[i][j] * initial_guess[j]) / A[i][xi]
                    print(-A[i][j], "*", initial_guess[j], "/", A[i][xi], "=", (-A[i][j] * initial_guess[j]) / A[i][xi])
                    iterations.arr[i] = round(iterations.arr[i], 14)
            xi = (xi + 1) % matrix.row
        pprint(iterations.arr)
        initial_guess = iterations.arr
        iterations.arr = [0] * matrix.row

def gauss_seidel(matrix, initial_guess, tol):
    return
    
c = ConvolutionalMatrix(n = 5)
x_stream = c.gen_random_x_stream()
y_stream = c.gen_y_stream(x_stream)
pprint(y_stream.arr)

m = Matrix(3, 3)
m.mat = [[101, 21, -1], [11, 8, 3], [-2, -1, 10]]
m.b = [7, -4, 3]
v = Vector(3)
v.arr = [0, 0, 0]
jacobi(m, [0, 0, 0], 5)
