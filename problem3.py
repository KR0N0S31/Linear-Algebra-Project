class Matrix:
    def __init__(self, n, m):
        self.row = n
        self.col = m
        self.mat = [[0 for y in range(self.col)] for x in range(self.row)]

    def multiply_by_vector(self, vector):
        result = Vector(self.row)
        for i in range(self.row):
            for j in range(self.col):
                result.arr[i] += self.mat[i][j] * vector.arr[j]
        return result

class Vector:
    def __init__(self, n):
        self.row = n
        self.col = 1
        self.arr = [0] * n

    def fill(self, entry):
        for i in range(self.row):
            self.arr[i] = entry

    def equals_zero_vector(self):
        zero_count = 0
        for i in range(self.row):
            if self.arr[i] == 0:
                zero_count += 1
        return zero_count == self.row

def power_method(A, u0, tol):
    if u0.equals_zero_vector():
        print("Error: u0 must be a non-zero vector.")
        return
    w = Vector(A.row)
    w.arr[0] = 1
    for i in range(5):
        v = A.multiply_by_vector(u0)
        print(v.arr)
        print(v.arr[0] / u0.arr[0])
        u0 = v

A = Matrix(2, 2)
A.mat = [[3, 1], [2, 4]]
u0 = Vector(A.row)
u0.arr = [1, 21]
power_method(A, u0, 0)
