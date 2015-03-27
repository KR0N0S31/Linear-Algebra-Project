import random
import copy
from pprint import pprint

class Matrix:
    def __init__(self, n, m):
        self.row = n
        self.col = m
        self.mat = [[0 for y in range(self.col)] for x in range(self.row)]
        self.b   = Vector(self.row)

    def augment_with(self, b):
        self.b = b

    def crop(self, n, m):
        result = Matrix(n, m)
        for i in range(n):
            for j in range(m):
                result.mat[i][j] = self.mat[i][j]
        return result

class Vector:
    def __init__(self, n):
        self.row = n
        self.col = 1
        self.arr = [0] * n

    def print_as_bits(self):
        print("[", end="")
        for i in range(self.row):
            if i == self.row - 1:
                print("%d" % self.arr[i], end="")
            else:
                print("%d," % self.arr[i], end=" ")
        print("]")

class ConvolutionalMatrix:
    def __init__(self, n):
        self.n = n
        self.row = self.n + 3
        self.col = self.n
        self.A0 = Matrix(self.row, self.col)
        self.A1 = Matrix(self.row, self.col)

    def gen_random_x_stream(self):
        binary_stream = Vector(self.n)
        random.seed()
        for i in range(self.n):
            binary_stream.arr[i] = random.randint(0, 1)
        return binary_stream

    def gen_y_stream(self, x_stream):
        jth_pos = 0
        # Generate A0 and A1 simultaneously
        for x in range(self.row):
            for y in range(self.col):
                # Formula for generating A0
                if y == jth_pos or y == jth_pos - 2 or y == jth_pos - 3:
                    self.A0.mat[x][y] = 1
                else:
                    self.A0.mat[x][y] = 0
                # Formula for generating A1
                if y == jth_pos or y == jth_pos - 1 or y == jth_pos - 3:
                    self.A1.mat[x][y] = 1
                else:
                    self.A1.mat[x][y] = 0
            jth_pos += 1
        y0 = self.multiply_binary_matrix_by_vector(self.A0, x_stream)
        y1 = self.multiply_binary_matrix_by_vector(self.A1, x_stream)
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

    def decode_y_stream(self, y_stream):
        # The x byte stream can be derived from either A0|y0 or A1|y1
        # A0 needs to be an n x (n + 1) matrix so it needs to be "cropped"
        # y0 also needs to be "cropped"
        cropped_y0 = Vector(self.col)
        for i in range(self.col):
            y_str = y_stream.arr[i]
            y0_num = y_str[0]
            cropped_y0.arr[i] = float(y0_num)
        cropped_A0 = self.A0.crop(self.col, self.col)
        cropped_A0.augment_with(cropped_y0)
        v = Vector(self.col)
        gauss_seidel(cropped_A0, v, 0.000001, 1)
        jacobi(cropped_A0, v, 0.000001, 1)
        
def jacobi(matrix, current_guess, tol, decode_binary_stream):
    A = matrix.mat
    b = matrix.b
    new_guess = Vector(matrix.row)
    xi = 0
    num_iteration = 0
    current_error_tol = float("inf")
    while (current_error_tol > tol and num_iteration < 100):
        for i in range(matrix.row):
            for j in range(matrix.col):
                if j == xi:
                    new_guess.arr[i] += b.arr[i] / A[i][xi]
                    new_guess.arr[i] = round(new_guess.arr[i], 14)
                    if decode_binary_stream:
                        new_guess.arr[i] %= 2
                else:
                    new_guess.arr[i] += (-A[i][j] * current_guess.arr[j]) / A[i][xi]
                    new_guess.arr[i] = round(new_guess.arr[i], 14)
                    if decode_binary_stream:
                        new_guess.arr[i] %= 2
            xi = (xi + 1) % matrix.col
        current_error_tol = get_error_tolerance(new_guess, current_guess)
        current_guess.arr = new_guess.arr
        new_guess.arr = [0] * matrix.row
        num_iteration += 1
    if num_iteration == 100:
        print("Failed to converge with Jacobi method after 100 iterations. Current solution:", current_guess.arr)
        return num_iteration
    if decode_binary_stream:
        print("Jacobi solution:", end=" ")
        current_guess.print_as_bits()
    else:
        print("Jacobi solution:", current_guess.arr)
    print("Iterations using Jacobi:", num_iteration)
    return num_iteration
        
def gauss_seidel(matrix, current_guess, tol, decode_binary_stream):
    A = matrix.mat
    b = matrix.b
    new_guess = Vector(matrix.row)
    xi = 0
    num_iteration = 0
    current_error_tol = float("inf")
    while (current_error_tol > tol and num_iteration < 100):
        # Current iteration before changes to accurately compute error
        _current_guess = copy.deepcopy(current_guess) 
        for i in range(matrix.row):
            for j in range(matrix.col):
                if j == xi:
                    new_guess.arr[i] += b.arr[i] / A[i][xi]
                    new_guess.arr[i] = round(new_guess.arr[i], 14)
                    if decode_binary_stream:
                        new_guess.arr[i] %= 2
                else:
                    new_guess.arr[i] += (-A[i][j] * current_guess.arr[j]) / A[i][xi]
                    new_guess.arr[i] = round(new_guess.arr[i], 14)
                    if decode_binary_stream:
                        new_guess.arr[i] %= 2
                current_guess.arr[i] = new_guess.arr[i]
            xi = (xi + 1) % matrix.col
        current_error_tol = get_error_tolerance(new_guess, _current_guess)
        current_guess.arr = new_guess.arr
        new_guess.arr = [0] * matrix.row
        num_iteration += 1
    if num_iteration == 100:
        print("Failed to converge with Jacobi method after 100 iterations. Current solution:", current_guess.arr)
        return num_iteration
    if decode_binary_stream:
        print("Gauss-Seidel solution:", end=" ")
        current_guess.print_as_bits()
    else:
        print("Gauss-Seidel solution:", current_guess.arr)
    print("Iterations using Gauss-Seidel:", num_iteration)
    return num_iteration
        

def get_error_tolerance(new, old):
    # Error defined as ||x^n - x^(n + 1)||
    error = Vector(new.row)
    for i in range(new.row):
        error.arr[i] = abs(old.arr[i] - new.arr[i])
    return max(error.arr)
    
c = ConvolutionalMatrix(n = 10)
x_stream = c.gen_random_x_stream()
print("x_stream =", x_stream.arr)
#x_stream.arr = [1, 0, 1, 1, 0]
y_stream = c.gen_y_stream(x_stream)
c.decode_y_stream(y_stream)

m = Matrix(5, 5)
m.mat = [[1, 0, 0, 0, 0], [0, 1, 0, 0, 0], [1, 0, 1, 0, 0], [1, 1, 0, 1, 0], [0, 1, 1, 0, 1]]
b = Vector(4)
b.arr = [1, 0, 0, 0, 1]
m.augment_with(b)

v = Vector(5)
v.arr = [0, 0, 0, 0, 0]
#jacobi(m, v, 0.01, 1)
#gauss_seidel(m, v, 0.01, 1)
