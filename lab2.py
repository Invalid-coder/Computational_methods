import math
import sympy
from math import log2
import matplotlib.pyplot as plt
from sympy import diff, symbols
from sympy.solvers import solve
from sympy.parsing.sympy_parser import parse_expr
import numpy as np

n = [1, 2, 5, 10, 20, 33, 34, 35]#, 100]


def f(x):
    return log2(3 + x)


def fact(n):
    res = 1

    for i in range(2, n + 1):
        res *= i

    return res


def Qnk(k, x, nodes):
    res = 1

    for i in range(k):
        res *= (x - nodes[i])

    return res


def Rnk(k, n, x, nodes):
    res = 1

    for i in range(n + 1):
        if i == k:
            continue

        res *= (nodes[k] - nodes[i])

    return res


def funk(k, x, nodes):
    res = 0

    for i in range(k + 1):
        res += f(nodes[i]) / Rnk(i, k, x, nodes)

    return res


def Pn(x, nodes):
    res = 0

    for k in range(len(nodes)):
        res += funk(k, x, nodes) * Qnk(k, x, nodes)

    return res


def U(x, nodes):
    res = 1

    for i in range(len(nodes)):
        res *= (x - nodes[i])

    return res


def main(points, nodes):
    for i in range(len(nodes)):
        y = [Pn(p, nodes[i]) for p in points]
        print(y)
        plt.plot(points, y, label='n = {}'.format(n[i]))
    plt.xlabel("x")
    plt.ylabel("Pn(x)")
    plt.title("Polinom")
    plt.legend()
    plt.show()


def errors(points, nodes):
    x = symbols('x')
    f = parse_expr("log2(3 + x)", local_dict={"log2": lambda x: sympy.log(x, 2)})
    for i in range(len(nodes)):
        d = diff(f, x, n[i] + 1)
        y = [(d.subs(x, p) * U(p, nodes[i])) / fact(n[i] + 1) for p in points]
        plt.plot(points, y, label='n = {}'.format(n[i]))
    plt.title("Errors")
    plt.legend()
    plt.show()


def print_table(nodes):
    prev = [f(x) for x in nodes]
    n = len(nodes)
    print(prev)
    for i in range(1, n):
        current = [prev[j] - prev[j - 1] for j in range(1, len(prev))]
        prev = current
        print(current)

def analitically(n, a, b, nodes):
    x = symbols('x')
    f = parse_expr("log2(3 + x)", local_dict={"log2": lambda x: sympy.log(x, 2)})
    f = diff(f, x, n + 1)

    critical_points = solve(diff(f, x), x)
    critical_points.extend([a, b])

    M_max = max(map(lambda val: abs(f.subs(x, val)), critical_points))
    U = 1

    for i in range(len(nodes)):
        U *= (x - nodes[i])

    critical_points = solve(diff(U, x), x)
    critical_points.extend([a, b])

    U_max = max(map(lambda val: abs(U.subs(x, val)), critical_points))
    alpha = (M_max * U_max) / fact(n + 1)
    return alpha


if __name__ == '__main__':
    a, b = -1, 1
    h = [(b - a) / n[i] for i in range(len(n))]
    points = np.linspace(a, b, num=100)
    nodes = []

    for i in range(len(n)):
        row = [a]
        for j in range(1, n[i] + 1):
            row.append(row[j - 1] + h[i])
        nodes.append(row)

    main(points, nodes)
    errors(points, nodes)
    print_table(nodes[3])
    print()
    print(analitically(n[3], a, b, nodes[3]))

    nodes = []

    for i in range(len(n)):
        row = [a]
        for m in range(1, n[i] + 1):
            row.append(math.cos(((2*m - 1) * math.pi) / (2 * n[i])))
        nodes.append(row)

    main(points, nodes)
    errors(points, nodes)
    print_table(nodes[3])
    print()
    print(analitically(n[3], a, b, nodes[3]))