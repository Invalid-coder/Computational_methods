import sympy
from math import log2
from math import log
from sympy import diff, symbols
from sympy.solvers import solve
from sympy.parsing.sympy_parser import parse_expr

n = [n for n in range(1, 11)]

X = {1: [0.5],
     2: (0.21132487, 0.78867513),
     3: (0.11270167, 0.5, 0.88729833),
     4: (0.06943184, 0.33000948, 0.66999052, 0.93056816),
     5: (0.04691008, 0.23076534, 0.5, 0.76923466, 0.95308992)
    }

consts = {1: [1],
          2: (0.5, 0.5),
          3: (5/18, 4/9, 5/18),
          4: (0.17392742, 0.32607258, 0.32607258, 0.17392742),
          5: (0.11846344, 0.23931433, 0.28444444, 0.23931433, 0.11846344)
         }


def f(x):
    return log2(3 + x)


def trapeze_formula(n, x, a, b):
    res = f(x[0])

    for i in range(1, n - 1):
        res += 2 * f(x[i])

    res += f(x[n - 1])

    return (b - a) * res / (2 * n)


def trapeze_error(n, a, b):
    x = symbols('x')
    f = parse_expr("log2(3 + x)", local_dict={"log2": lambda x: sympy.log(x, 2)})
    f = diff(f, x, 2)

    critical_points = solve(diff(f, x), x)
    critical_points.extend([a, b])

    M_max = max(map(lambda val: abs(f.subs(x, val)), critical_points))

    alpha = ((b - a) ** 3) * M_max / (12 * (n ** 2))
    return alpha


def simpson_formula(n, x, a, b):
    res = f(x[0])

    for i in range(1, n - 1):
        if i % 2 == 1:
            res += 4 * f(x[i])
        else:
            res += 2 * f(x[i])

    res += f(x[n - 1])

    return (b - a) * res / (6 * n)


def simpson_error(n, a, b):
    x = symbols('x')
    f = parse_expr("log2(3 + x)", local_dict={"log2": lambda x: sympy.log(x, 2)})
    f = diff(f, x, 4)
    h = (b - a) / (2 * n)
    critical_points = solve(diff(f, x), x)
    critical_points.extend([a, b])

    M_max = max(map(lambda val: abs(f.subs(x, val)), critical_points))

    alpha = (b - a) * (h ** 4) * M_max / 180
    return alpha


def gauss_formula(n, a, b):
    assert n <= 5

    x = X[n]
    C = consts[n]
    res = 0

    for i in range(n):
        z = a + (b - a) * x[i]
        res += C[i] * f(z)

    return (b - a) * res


def precise_value(n, a, b):
    F_b = (3 + b) * log(3 + b) / log(2) - (3 + b) / log(2)
    F_a = (3 + a) / log(2) - (3 + a) * log(3 + a) / log(2)

    return F_b - F_a


def create_table(nodes, a, b):
    for i, x in enumerate(nodes):
        print("n = ", n[i])
        print("Nodes: ", x)
        print("Precise value = ", precise_value(n[i], a, b))
        print("Trapeze formula = ", trapeze_formula(n[i], nodes[i], a, b))
        print("Simpson formula = ", simpson_formula(n[i], nodes[i], a, b))

        if n[i] <= 5:
            print("Gauss formula = ", gauss_formula(n[i], a, b))

        print("Trapeze error = ", trapeze_error(n[i], a, b))
        print("Simposon error = ", simpson_error(n[i], a, b))


if __name__ == '__main__':
    a, b = -1, 1
    h = [(b - a) / n[i] for i in range(len(n))]
    nodes = []

    for i in range(len(n)):
        row = [a]
        for j in range(1, n[i] + 1):
            row.append(row[j - 1] + h[i])
        nodes.append(row)

    create_table(nodes, a, b)