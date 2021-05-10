import math


def a(i, j):
    return math.log(i) / (i + j + 1)


def b(i):
    return 0.1 * math.exp(-i)


def get_matrix(n, m):
    A = []
    for i in range(1, n + 1):
        A.append([])
        for j in range(1, m + 1):
            A[i - 1].append(a(i, j))
        A[i - 1].append(b(i))
    return A


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


def dot_product(a, b):
    return sum(a[i] * b[i] for i in range(len(a)))


def length_sqr(v):
    return sum(v[i] * v[i] for i in range(len(v)))


def orthogonal_vectors(A):
    r, f = [], []
    n = len(A)
    for i in range(n):
        r.append(A[i][0:-1])
        f.append(A[i][-1])
    R, F = [r[0]], [f[0]]
    for i in range(1, n):
        res = r[i]
        for j in range(i):
            res = subtract(res, multiply(R[j], dot_product(r[i], R[j]) / length_sqr(R[j])))
        R.append(res)
    for i in range(1, n):
        res = f[i]
        for j in range(i):
            res -= F[j] * dot_product(r[i], R[j]) / length_sqr(R[j])
        F.append(res)
    res = multiply(R[0], F[0] / length_sqr(R[0]))
    for i in range(1, n):
        res = add(res, multiply(R[i], F[i] / length_sqr(R[i])))
    return res


if __name__ == '__main__':
    n = 5
    A = get_matrix(n, n + 1)#[[2, 3, 4, 20], [3, 2, 1, 10], [4, 1, 2, 12]] homework matrix
    #A1 = [[2, 3, -4, 1], [-3, -1, 5, 1], [4, 3, -6, 1]] #classwork matrix
    for x in A:
        print(x)
    print(gaussian_method(A))
    print(orthogonal_vectors(A))