import math


def get_matrix(n):
    A = [[0] * n for _ in range(n)]
    b = [[1 / (i * i)] for i in range(1, n + 1)]
    for i in range(1, n + 1):
        A[i - 1][i - 1] = 2
        if i != n: A[i - 1][i] = 1 / (i * math.log(i + 1))
    return A, b


def get_B_f(A, b):
    B = [[-A[i][j] / A[i][i] if i != j else 0 for j in range(len(A[0]))] for i in range(len(A))]
    f = [[b[i][0] / A[i][i]] for i in range(len(b))]
    return B, f


def convergence_check(B):
    return max(list(sum(map(lambda x: abs(x), B[i])) for i in range(len(B)))) < 1


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


def add_vectors(a, b):
    assert len(a) == len(b)
    return [[a[i][0] + b[i][0]] for i in range(len(a))]


def subtract_vectors(a, b):
    assert len(a) == len(b)
    return [[a[i][0] - b[i][0]] for i in range(len(a))]


def multiply(a, v):
    return [[a * x[0]] for x in v]


def length(A, b, X):
    r = []
    for i in range(len(A)):
        r_i = - b[i][0]
        for j in range(len(A[0])):
            r_i += A[i][j] * X[j][0]
        r.append(r_i)
    return math.sqrt(sum(map(lambda x: x * x, r)))


def simple_iteration_method(A, b, eps):
    B, f = get_B_f(A, b)
    X = f.copy()
    i = 0
    if convergence_check(B):
        print("Method of simple iteration is convergent")
        while length(A, b, X) > eps:
            X = add_vectors(matrix_multiplication(B, X), f)
            i += 1
    else:
        print("Method of simple iteration is not convergent")
        return None
    return X, i, length(A, b, X)


def seidel_method(A, b, eps):
    B, f = get_B_f(A, b)
    X = [[f[i][0]] for i in range(len(f))]
    i = 0
    if convergence_check(B):
        print("Seidel's method is convergent")
        while length(A, b, X) > eps:
            for i in range(len(B)):
                x_i = f[i][0]
                for j in range(len(B[0])):
                    x_i += B[i][j] * X[j][0]
                X[i][0] = x_i
            i += 1
    else:
        print("Seidel's method is not convergent")
        return None
    return X, i, length(A, b, X)


def convergence_check_gd(A):
    min_diag = min(A[i][i] for i in range(len(A)))
    return all([A[i][j] < min_diag for i in range(len(A)) for j in range(len(A[0])) if i != j])


def transposed(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def dot_product(a, b):
    return sum(a[i][0] * b[i][0] for i in range(len(a)))


def vector_length(v):
    return math.sqrt(sum(x[0] * x[0] for x in v))


def gradient_descent(A, b, eps):
    A_t = transposed(A)
    X = [[b[i][0] / A[i][i]] for i in range(len(A))]
    i = 0
    if convergence_check_gd(A):
        print("Gradient_descent method is convergent")
        r = subtract_vectors(matrix_multiplication(A, X), b)
        while vector_length(r) > eps:
            product = matrix_multiplication(matrix_multiplication(A, A_t), r)
            delta = dot_product(r, product) / dot_product(product, product)
            X = subtract_vectors(X, multiply(delta, matrix_multiplication(A_t, r)))
            r = subtract_vectors(matrix_multiplication(A, X), b)
            i += 1
    else:
        print("Gradient_descent is not convergent")
        return None
    return X, i, vector_length(r)


if __name__ == '__main__':
    n = 100
    eps = 10**(-5) #0.4
    A, b = get_matrix(n)
    #A = [[8, 3, 2], [2, 10, 4], [5, 2, 8]]
    #b = [[2], [1], [2]]
    X1, i1, r1 = simple_iteration_method(A, b, eps)
    print("X = ", X1)
    print("interations_number = ", i1)
    print("precison = ", r1)
    X2, i2, r2 = seidel_method(A, b, eps)
    print("X = ", X2)
    print("interations_number = ", i2)
    print("precison = ", r2)
    X3, i3, r3 = gradient_descent(A, b, eps)
    print("X = ", X3)
    print("interations_number = ", i3)
    print("precison = ", r3)



