"""
https://en.wikipedia.org/wiki/Divisor_function
"""
from math import exp, log
from divisor import divisors
import matplotlib.pyplot as plt
import matplotlib.style as style


def harmonic(n: int) -> float:
    """Computes the nth harmonic number"""
    h = 0.0
    for i in range(1, n+1):
        h += 1.0/i
    return h


def robin_ineq(n):
    lhs = divisors(1, n)
    h_n = harmonic(n)
    rhs = h_n + exp(h_n) * log(h_n)
    return rhs - lhs


def plotting(ls_x, ls_y):
    style.use('ggplot')  # Set the plot style

    fig, ax = plt.subplots(figsize=(8, 6))  # Set the plot size

    ax.scatter(ls_x, ls_y, s=2, label='$h_n + e^{h_n}\ln(h_n) - \sigma(n)$')  # Use LaTeX to render label

    ax.axvline(x=5040, ls='--', c='red', linewidth=1, label='starting value = 5040')  # Add vertical line at n=5040

    ax.axhline(y=0, color='black', linestyle='--')  # Add horizontal line at y=0

    ax.set_title("Robin's inequality", fontsize=20, fontweight='bold')  # Set the plot title
    ax.set_xlabel("$n$", fontsize=16)  # Use LaTeX to render axis label
    ax.set_ylabel("$h_n + e^{h_n}\ln(h_n) - \sigma(n)$", fontsize=16)  # Use LaTeX to render axis label
    ax.tick_params(axis='both', which='major', labelsize=12)  # Set the font size of the tick labels

    ax.legend(fontsize=12)  # Add the legend
    plt.tight_layout()  # Make sure the plot is well-aligned

    plt.show()


if __name__ == '__main__':
    n_range = range(5041, 10000)  # change upper bound to 20,000 to achieve the stored plot
    ls = [robin_ineq(n) for n in n_range]

    plotting(ls_x=n_range, ls_y=ls)
