import math
import copy
from pprint import pprint

class Matrix:
    def __init__(self, A = None, b = None):
        self.row = len(A) if A != None else 0
        self.col = len(A[0]) if A != None else 0
        self.__A = A

    def set_A(self, matrix):
        self.__A = matrix
        self.row = len(self.__A)
        self.col = len(self.__A[0])

    def get_A(self):
        return self.__A

    def multiply_by_vector(self, vector):
        result = [0] * self.row
        for i in range(self.row):
            for j in range(self.col):
                result[i] += self.__A[i][j] * vector[j]
        return result

def power_method(A, u0, tol):
    if equals_zero_vector(u0):
        print("Error: u0 must be a non-zero vector.")
        return
    elif A.row != A.col:
        print("Error: A must be of size n x n")
        return
    curr_eigen_val = 0
    prev_eigen_val = 0
    curr_eigen_vec = [0] * len(u0)
    w = [0] * A.row
    w[0] = 1
    k = float("-inf")
    i = 0
    while k < tol:
        v = A.multiply_by_vector(u0)
        prev_eigen_val = curr_eigen_val
        curr_eigen_val = v[0] / u0[0]
        curr_eigen_vec = divide_vec_by_norm(u0)
        if i > 1:
            r = abs(curr_eigen_val / prev_eigen_val) # convergence factor
            k = math.log1p(tol) / math.log1p(r) # current tolerance
        u0 = v
        i += 1
    print("After %d iterations, the following values were found:" % i)
    print("Eigen value =", curr_eigen_val, "\nNormalized eigen vectors =", curr_eigen_vec)
    print()
    return curr_eigen_val, curr_eigen_vec, i

def equals_zero_vector(vector):
        zero_count = 0
        for i in range(len(vector)):
            if vector[i] == 0:
                zero_count += 1
        return zero_count == len(vector)

def divide_vec_by_norm(vector):
    norm = 0
    for i in range(len(vector)):
        norm += vector[i] * vector[i]
    norm = math.sqrt(norm)
    result = [0] * len(vector)
    for i in range(len(vector)):
        result[i] = vector[i] / norm
    return result

def compute_leslie_populations(Leslie, x0):
    x1 = [0] * Leslie.row
    distributions = [[0 for x in range(Leslie.row)] for x in range(6)]
    distributions[0] = x0
    for k in range(5):
        for i in range(Leslie.row):
            for j in range(Leslie.col):
                x1[i] += Leslie.get_A()[i][j] * x0[j]
            x1[i] = round(x1[i], 0)
            distributions[k + 1][i] = x1[i] 
        x0 = x1
    prev_pop = 0
    total_pop = 0
    for i in range(len(distributions)):
        for j in range(len(distributions[0])):
            total_pop += distributions[i][j]
        year = 2000 + (i * 10)
        print("Population distribution for %d =" % year)
        pprint(distributions[i])
        print("Total population for %d = %d" % (year, total_pop))
        if i > 0:
            percent_change = ((total_pop - prev_pop) / total_pop) * 100
            print("Percent change in population from", year - 10, "=", percent_change, "%")
        print()
        prev_pop = total_pop
        total_pop = 0  

Leslie = Matrix([[0, 1.2, 1.1, .9, .1, 0, 0, 0, 0],
                 [0.7, 0, 0, 0, 0, 0, 0, 0, 0],
                 [0, .85, 0, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0.9, 0, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0.9, 0, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0.88, 0, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0.8, 0, 0, 0],
                 [0, 0, 0, 0, 0, 0, 0.77, 0, 0],
                 [0, 0, 0.9, 0, 0, 0, 0, 0.4, 0]])
x0 = [2.1, 1.9, 1.8, 2.1, 2.0, 1.7, 1.2, 0.9, 0.5]
for i in range(len(x0)):
    x0[i] *= 100000
    
# Problem 3b    
compute_leslie_populations(Leslie, x0)
# Problem 3c
power_method(A = Leslie, u0 = [1, 1, 1, 1, 1, 1, 1, 1, 1], tol = 1)

# Problem 3d
Leslie.get_A()[0][1] /= 2
compute_leslie_populations(Leslie, x0)
power_method(A = Leslie, u0 = [1, 1, 1, 1, 1, 1, 1, 1, 1], tol = 1)
