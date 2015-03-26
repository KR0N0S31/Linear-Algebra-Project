import random
from pprint import pprint

class Matrix:
    def __init__(self, n, m):
        self.row = n
        self.col = m
        self.arr = [[0 for y in range(self.col)] for x in range(self.row)]

    def augment_with(self, b):
        new_arr = [[0 for y in range(self.col + 1)] for x in range(self.row)]
        for x in range(self.row):
            for y in range(self.col + 1):
                if y == self.col:
                    new_arr[x][y] = b.arr[x]
                else:
                    new_arr[x][y] = self.arr[x][y]
        self.col += 1
        self.arr = new_arr

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
                    A0.arr[x][y] = 1
                else:
                    A0.arr[x][y] = 0
                # Formula for generating A1
                if y == jth_pos or y == jth_pos - 1 or y == jth_pos - 3:
                    A1.arr[x][y] = 1
                else:
                    A1.arr[x][y] = 0
            jth_pos += 1
        y0 = self.multiply_binary_matrix_by_vector(A0, x_stream)
        y1 = self.multiply_binary_matrix_by_vector(A1, x_stream)
        y  = self.merge_y0_with_y1(y0, y1)
        return y

    def multiply_binary_matrix_by_vector(self, matrix, vector):
        result = Vector(matrix.row)
        for i in range(matrix.row):
            for j in range(vector.row):
                result.arr[i] += matrix.arr[i][j] * vector.arr[j]
                result.arr[i] %= 2
        return result

    def merge_y0_with_y1(self, y0, y1):
        y = Vector(y0.row)
        for x in range(y0.row):
            concat = str(y0.arr[x]) + str(y1.arr[x])
            y.arr[x] = concat
        return y

def jacobi(matrix, initial_guess, tol):
    return

def gauss_seidel(matrix, initial_guess, tol):
    return
    
c = ConvolutionalMatrix(n = 5)
x_stream = c.gen_random_x_stream()
y_stream = c.gen_y_stream(x_stream)
pprint(y_stream.arr)

