import random
import copy
from pprint import pprint

class Matrix:
    def __init__(self, mat = None, b = None):
        self.row = len(mat) if mat != None else 0
        self.col = len(mat[0]) if mat != None else 0
        self.__mat = mat
        self.__b = b

    def set_matrix(self, matrix):
        self.__mat = matrix
        self.row = len(self.__mat)
        self.col = len(self.__mat[0])

    def set_b(self, b):
        self.__b = b

    def get_mat(self):
        return self.__mat

    def get_b(self):
        return self.__b

    def crop(self, n, m):
        result = Matrix(mat = [[0 for x in range(m)] for x in range(n)])
        for i in range(n):
            for j in range(m):
                result.__mat[i][j] = self.__mat[i][j]
        return result

class ConvolutionalMatrix:
    def __init__(self, n):
        self.n = n
        self.row = self.n + 3
        self.col = self.n
        m = [[0 for x in range(self.col)] for x in range(self.row)]
        self.A0 = Matrix(mat = m)
        self.A1 = Matrix(mat = m)

    def gen_random_x_stream(self):
        binary_stream = [0] * self.n
        random.seed()
        for i in range(self.n):
            binary_stream[i] = random.randint(0, 1)
        print("randomly generated x_stream of length", self.n, "=", binary_stream)
        return binary_stream

    def gen_y_stream(self, x_stream):
        jth_pos = 0
        # Generate A0 and A1 simultaneously
        for x in range(self.row):
            for y in range(self.col):
                # Formula for generating A0
                if y == jth_pos or y == jth_pos - 2 or y == jth_pos - 3:
                    self.A0.get_mat()[x][y] = 1
                else:
                    self.A0.get_mat()[x][y] = 0
                # Formula for generating A1
                if y == jth_pos or y == jth_pos - 1 or y == jth_pos - 3:
                    self.A1.get_mat()[x][y] = 1
                else:
                    self.A1.get_mat()[x][y] = 0
            jth_pos += 1
        y0 = self.multiply_binary_matrix_by_vector(self.A0, x_stream)
        y1 = self.multiply_binary_matrix_by_vector(self.A1, x_stream)
        y  = self.merge_y0_with_y1(y0, y1)
        print("A0 =")
        pprint(self.A0.get_mat())
        print("A1 =")
        pprint(self.A1.get_mat())
        print("y0 =", y0)
        print("y1 =", y1)
        print("y_stream =", y)
        print()
        return y

    def multiply_binary_matrix_by_vector(self, matrix, vector):
        result = [0] * matrix.row
        for i in range(matrix.row):
            for j in range(matrix.col):
                result[i] += matrix.get_mat()[i][j] * vector[j]
                result[i] %= 2
        return result

    def merge_y0_with_y1(self, y0, y1):
        y = [0] * len(y0)
        for x in range(len(y0)):
            concat = str(y0[x]) + str(y1[x])
            y[x] = concat
        return y

    def decode_y_stream(self, y_stream):
        # The x byte stream can be derived from either A0|y0 or A1|y1
        # A0 needs to be an n x (n + 1) matrix so it needs to be "cropped"
        # y0 also needs to be "cropped"
        cropped_y0 = [0] * self.col
        for i in range(self.col):
            y_str = y_stream[i]
            y0_num = y_str[0]
            cropped_y0[i] = float(y0_num)
        cropped_A0 = self.A0.crop(self.col, self.col)
        cropped_A0.set_b(cropped_y0)
        v = [0] * self.col
        gauss_seidel(cropped_A0, v, 0.000001, 1)
        jacobi(cropped_A0, v, 0.000001, 1)

# This method makes two assumptions:
#   (1) that the system given by a matrix  has a unique solution and
#   (2) that the coefficient matrix A has no zeros on its main diagonal
def jacobi(matrix, guess, tol, decode_binary_stream):
    A = matrix.get_mat()
    b = matrix.get_b()
    new_guess = [0] * matrix.row
    xi = 0
    num_iteration = 0
    current_error_tol = float("inf")
    while (current_error_tol > tol and num_iteration < 100):
        for i in range(matrix.row):
            for j in range(matrix.col):
                if j == xi:
                    if A[i][xi] == 0:
                        new_guess[i] += 0
                    else:
                        new_guess[i] += b[i] / A[i][xi]
                    new_guess[i] = round(new_guess[i], 14)
                    if decode_binary_stream:
                        new_guess[i] %= 2
                else:
                    if A[i][xi] == 0:
                        new_guess[i] += 0
                    else:
                        new_guess[i] += (-A[i][j] * guess[j]) / A[i][xi]
                    new_guess[i] = round(new_guess[i], 14)
                    if decode_binary_stream:
                        new_guess[i] %= 2
            xi = (xi + 1) % matrix.col
        current_error_tol = get_error_tolerance(new_guess, guess)
        guess = new_guess
        new_guess = [0] * matrix.row
        num_iteration += 1
    if num_iteration == 100:
        print("Failed to converge with Jacobi method after 100 iterations.\nCurrent solution:", guess.__v)
        return num_iteration
    if decode_binary_stream:
        print("Jacobi solution:", end=" ")
        print_as_bits(guess)
    else:
        print("Jacobi solution:", guess)
    print("Iterations using Jacobi:", num_iteration)
    print()
    return num_iteration
        
def gauss_seidel(matrix, guess, tol, decode_binary_stream):
    A = matrix.get_mat()
    b = matrix.get_b()
    new_guess = [0] * matrix.row
    xi = 0
    num_iteration = 0
    current_error_tol = float("inf")
    while (current_error_tol > tol and num_iteration < 100):
        # Current iteration before changes to accurately compute error
        _guess = copy.deepcopy(guess) 
        for i in range(matrix.row):
            for j in range(matrix.col):
                if j == xi:
                    new_guess[i] += b[i] / A[i][xi]
                    new_guess[i] = round(new_guess[i], 14)
                    if decode_binary_stream:
                        new_guess[i] %= 2
                else:
                    new_guess[i] += (-A[i][j] * guess[j]) / A[i][xi]
                    new_guess[i] = round(new_guess[i], 14)
                    if decode_binary_stream:
                        new_guess[i] %= 2
                guess[i] = new_guess[i]
            xi = (xi + 1) % matrix.col
        current_error_tol = get_error_tolerance(new_guess, _guess)
        guess = new_guess
        new_guess = [0] * matrix.row
        num_iteration += 1
    if num_iteration == 100:
        print("Failed to converge with Gauss-Seidel method after 100 iterations.\nCurrent solution:", guess.__v)
        return num_iteration
    if decode_binary_stream:
        print("Gauss-Seidel solution:", end=" ")
        print_as_bits(guess)
    else:
        print("Gauss-Seidel solution:", guess)
    print("Iterations using Gauss-Seidel:", num_iteration)
    print()
    return num_iteration
        
def get_error_tolerance(new, old):
    # Error defined as ||x^n - x^(n + 1)||
    error = [0] * len(new)
    for i in range(len(new)):
        error[i] = abs(old[i] - new[i])
    return max(error)

def print_as_bits(v):
    print("[", end="")
    for i in range(len(v)):
        if i == len(v) - 1:
            print("%d" % v[i], end="")
        else:
            print("%d," % v[i], end=" ")
    print("]")
  
c = ConvolutionalMatrix(n = 5)
x_stream = c.gen_random_x_stream()
y_stream = c.gen_y_stream(x_stream)
c.decode_y_stream(y_stream)

##m = Matrix()
##m.set_matrix([[5, -2, 3], [-3, 2, 1], [2, -1, -7]])
##m.set_b([-1, 2, 3])
##jacobi(matrix = m, guess = [0, 0, 0, 0, 0], tol = 0.0001, decode_binary_stream = 0)
##gauss_seidel(matrix = m, guess = [0, 0, 0, 0, 0], tol = 0.0001, decode_binary_stream = 0)
