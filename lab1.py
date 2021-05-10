import sympy
from math import log2
import matplotlib.pyplot as plt
from sympy import diff, symbols
from sympy.solvers import solve
from sympy.parsing.sympy_parser import parse_expr
import numpy as np

n = [1, 2, 5, 10, 20, 100]


def f(x):
    return log2(3 + x)


def fact(n):
    res = 1

    for i in range(2, n + 1):
        res *= i

    return res


def Qnk(k, x, nodes):
    res = 1

    for i in range(len(nodes)):
        if i == k:
            continue

        res *= (x - nodes[i]) / (nodes[k] - nodes[i])

    return res


def Pn(x, nodes):
    res = 0

    for k in range(len(nodes)):
        res += f(nodes[k]) * Qnk(k, x, nodes)

    return res


def U(x, nodes):
    res = 1

    for i in range(len(nodes)):
        res *= (x - nodes[i])

    return res


def main(points, nodes):
    for i in range(len(nodes)):
        y = [Pn(p, nodes[i]) for p in points]
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
    points = np.linspace(a, b, num=10)
    nodes = []

    for i in range(len(n)):
        row = [a]
        for j in range(1, n[i] + 1):
            row.append(row[j - 1] + h[i])
        nodes.append(row)

    main(points, nodes)
    errors(points, nodes)
    print(analitically(n[3], a, b, nodes[3]))
