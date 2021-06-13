import math
import numpy as np
import matplotlib.pyplot as plt


def f1(x):
    return math.exp(-x)


def f2(x):
    return x**3 + 2*x - 1


def f3(x):
    return x**3 + 2*x - 1 - math.exp(-x)


def f3_derivative(x):
    return 3*x**2 + 2 + math.exp(-x)


def polinom(x):
    return 0.093307 * (x ** 3) + 0.171087 * (x ** 2) + 0.189761 * x + 0.190053


def binary_search(a, b, eps):
    left = a
    right = b

    while right - left > eps:
        mid = (left + right) / 2
        if polinom(a) * polinom(mid) > 0:
            left = mid
        if polinom(a) * polinom(mid) < 0:
            right = mid
        if polinom(mid) == 0:
            return mid

    return left


def chord_method(a, b, eps):
    x_prev = a
    x_prev_ = b
    x_curr = x_prev - f3(x_prev) * (x_prev_ - x_prev) / (f3(x_prev_) - f3(x_prev))

    while abs(x_curr - x_prev) > eps:
        x_prev, x_curr = x_curr, x_prev - f3(x_prev) * (x_curr - x_prev) / (f3(x_curr) - f3(x_prev))

    return x_curr


def newton_method(a, b, eps):
    x_prev = a
    x_curr = x_prev - f3(x_prev) / f3_derivative(x_prev)

    while abs(x_curr - x_prev) > eps:
        x_prev, x_curr = x_curr, x_curr - f3(x_curr) / f3_derivative(x_curr)

    return x_curr


if __name__ == '__main__':
    a, b = 0, 2
    eps = 10 ** (-6)
    print(binary_search(a, b, eps))
    x = np.linspace(a, b, num=10)
    y1 = [f1(x_i) for x_i in x]
    y2 = [f2(x_i) for x_i in x]
    plt.plot(x, y1, label='f1(x)')
    plt.plot(x, y2, label='f2(x)')
    plt.title("Graphics")
    plt.legend()
    plt.show()
    print(chord_method(a, b, eps))
    print(newton_method(a, b, eps))
