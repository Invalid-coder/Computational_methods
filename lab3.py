from math import log2
from math import log

n = [10, 11]#[1, 2, 5, 10, 20, 100]


def f(x):
    return log2(3 + x)


def second_derivative(x):
    return (-1) / (((3 + x) ** 2) * log(2))


def third_derivative(x):
    return 2 / (((3 + x) ** 3) * log(2))


def get_item(i, j, array):
    if len(array[i]) <= j:
        return array[i][-1]
    else:
        return array[i][j]


def get_delta_array(x):
    prev = [f(x_i) for x_i in x]
    delta_y = []
    for i in range(1, len(prev)):
        current = [prev[j] - prev[j - 1] for j in range(1, len(prev))]
        delta_y.append(current)
        prev = current
    return delta_y


def second_derivative_numerically(delta_y, task, i, h):
    if task == 1:
        return (1 / (h**2)) * (get_item(1, i, delta_y) + (i - 1) * get_item(2, i, delta_y))
    elif task == 2:
        return (1 / (h**2)) * (get_item(1, i, delta_y) + (i - 1) * get_item(2, i, delta_y) + \
               + (6 * (i**2) - 18 * i + 11) * get_item(3, i, delta_y) / 12 + \
               + (10 * (i**3) - 60 * (i**2) + 105 * i - 50) * get_item(4, i, delta_y) / 60)


def third_derivative_numerically(delta_y, task, i, h):
    if task == 1:
        return (1 / (h**3)) * get_item(2, i, delta_y)
    elif task == 2:
        return (1 / (h**3)) * (get_item(2, i, delta_y) + (2 * i - 3) * get_item(3, i, delta_y) / 2 + \
               + (10 * (i**2) - 40 * i + 35) * get_item(4, i, delta_y) / 20)


def create_table(nodes):
    for i, x in enumerate(nodes):
        delta_y = get_delta_array(x)
        print("Nodes: ", x)
        print("Table for task 1")
        for j in range(len(x)):
            approx_second = second_derivative_numerically(delta_y, 1, j, h[i])
            approx_third = third_derivative_numerically(delta_y, 1, j, h[i])
            precise_second = second_derivative(x[j])
            precise_third = third_derivative(x[j])
            print("Second derivative in node: ", x[j])
            print("Precise value: ", precise_second)
            print("Approximate value: ", approx_second)
            print("Absolute error: ", abs(precise_second - approx_second))
            print("Third derivative in node: ", x[j])
            print("Precise value: ", precise_third)
            print("Approximate value: ", approx_third)
            print("Absolute error: ", abs(precise_third - approx_third))
            print()
        print()
        print("Table for task 2")
        for j in range(len(x)):
            approx_second = second_derivative_numerically(delta_y, 2, j, h[i])
            approx_third = third_derivative_numerically(delta_y, 2, j, h[i])
            precise_second = second_derivative(x[j])
            precise_third = third_derivative(x[j])
            print("Second derivative in node: ", x[j])
            print("Precise value: ", precise_second)
            print("Approximate value: ", approx_second)
            print("Absolute error: ", abs(precise_second - approx_second))
            print("Third derivative in node: ", x[j])
            print("Precise value: ", precise_third)
            print("Approximate value: ", approx_third)
            print("Absolute error: ", abs(precise_third - approx_third))
            print()


if __name__ == '__main__':
    a, b = -1, 1
    h = [(b - a) / n[i] for i in range(len(n))]
    nodes = []

    for i in range(len(n)):
        row = [a]
        for j in range(1, n[i] + 1):
            row.append(row[j - 1] + h[i])
        nodes.append(row)

    create_table(nodes)