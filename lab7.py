import math
from sympy import diff, symbols


def get_matrix(n):
    A = [[0] * n for _ in range(n)]
    for i in range(1, n + 1):
        for j in range(1, n):
            A[i - 1][i - 1] = math.log(i) / (i + j + 1)
        A[i - 1][n - 1] = 0.1 * math.exp(-i)
    return A


def get_matrix2(n):
    A = [[0] * n for _ in range(n)]
    for i in range(1, n + 1):
        A[i - 1][i - 1] = 2
        if i != n: A[i - 1][i] = 1 / (i * math.log(i + 1))
        else: A[i - 1][-1] = 1 / (i * i)
    return A


def matrix_multiplication(A, B):
    rows_a, cols_a = len(A), len(A[0])
    rows_b, cols_b = len(B), len(B[0])
    res_matrix = [[0 for j in range(cols_b)] for i in range(rows_a)]
    if cols_a == rows_b:
        for i in range(rows_a):
            for j in range(cols_b):
                res_matrix[i][j] = sum(A[i][k] * B[k][j] for k in range(cols_a))
    else:
        print("columns of the first matrix must be equal to the rows of the second matrix")
        return None
    return res_matrix


def div(row, a):
    return [x / a for x in row]


def multiply(row, a):
    return [a * x for x in row]


def subtract(a, b):
    return [a[i] - b[i] for i in range(len(a))]


def add(a, b):
    return [a[i] + b[i] for i in range(len(a))]


def gaussian_method(A):
    n = len(A)
    for i in range(n):
        if A[i][i] == 0:
            continue
        A[i] = div(A[i],A[i][i])
        for j in range(i + 1, n):
            A[j] = subtract(A[j], multiply(A[i], A[j][i]))
    for i in range(n - 1, -1, -1):
        for j in range(i - 1, -1, -1):
            A[j] = subtract(A[j], multiply(A[i], A[j][i]))
    return [A[i][-1] for i in range(n)]


def krylov_method(A):
    C = [[[0] for i in range(len(A))]]
    C[0][0][0] = 1
    for i in range(1, len(A) + 1):
        C.append(matrix_multiplication(A, C[i - 1]))
    B = [[C[i][j][0] for i in range(len(C))] for j in range(len(C[0]))]
    p = gaussian_method(B)
    n = len(A)
    lambda_ = symbols('λ')
    D = lambda_ ** n
    for i in range(n - 1, -1, -1):
        D -= p[i] * (lambda_ ** i)
    D *= (-1) ** n
    return D


def danilevsky_method(A):
    j = 2
    a = A[-1][-j]
    n = len(A)
    m = len(A[0])
    last_row = A[-1].copy()
    while a == 0 and j < m:
        j += 1
        a = A[-1][-j]
    if a:
        for i in range(n):
            A[i][-j], A[i][-2] = A[i][-2], A[i][-j]
        for i in range(n):
            A[i][-2] = A[i][-2] / a
        B = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        B[-2] = [-A[-1][i] / a if i != m - 2 else 1 / a for i in range(m)]
        B_inv = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        B_inv[-2] = last_row
        A = matrix_multiplication(A, B)
        A = matrix_multiplication(B_inv, A)
        B = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        a = A[-2][-2]
        for j in range(m):
            if j == m - 3:
                B[1][j] = 1 / a
            elif j == m - 2:
                B[1][j] = -A[-2][j - 1] / a
            else:
                B[1][j] = -A[-2][j] / a
        B_inv = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        B_inv[1] = A[-2].copy()
        A = matrix_multiplication(A, B)
        A = matrix_multiplication(B_inv, A)
        B = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        a = A[1][0]
        B[0] = [-A[1][j] / a if j != 0 else 1 / a for j in range(m)]
        B_inv = [[1 if i == j else 0 for j in range(m)] for i in range(n)]
        B_inv[1] = A[1].copy()
        A = matrix_multiplication(A, B)
        A = matrix_multiplication(B_inv, A)
    else:
        return None
    p = A[0]
    lambda_ = symbols('λ')
    D = lambda_ ** n
    for i in range(m - 1, -1, -1):
        D -= p[m - i - 1] * (lambda_ ** i)
    D *= (-1) ** n
    return D


def power_method(A, eps):
    z = [[1] for i in range(len(A))]
    prev = 1
    z = matrix_multiplication(A, z)
    curr = max(z, key=lambda x: x[0])[0]
    z = [[x[0] / curr] for x in z]
    while abs(curr - prev) > eps:
        z = matrix_multiplication(A, z)
        prev = curr
        curr = max(z, key=lambda x: x[0])[0]
        z = [[x[0] / curr] for x in z]
    return curr


def transposed(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def dot_product(a, b):
    return sum(a[i][0] * b[i][0] for i in range(len(a)))


def vector_length(v):
    return math.sqrt(sum(x[0] * x[0] for x in v))


if __name__ == '__main__':
    A = get_matrix(3)
    #A = [[-1, 1, 3, -2], [0, 2, 3, -4], [1, -5, 3, 2], [-4, 6, 3, 1]]
    print(krylov_method(A))
    #A = [[2, 1, 4, 1], [3, 3, 2, -2], [4, 2, -1, 3], [5, -1, 4, 2]]
    A = get_matrix2(5)
    print(danilevsky_method(A))
    eps = 10 ** (-5)
    A = get_matrix2(4)
    print(power_method(A, eps))

