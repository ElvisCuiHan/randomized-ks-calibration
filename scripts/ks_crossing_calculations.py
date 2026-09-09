#!/usr/bin/env python3
"""Exact KS crossing calculations.

This file evaluates the Smirnov-Birnbaum-Tingey one-sided tail probability
using the log-sum-exp form described in the manuscript, and computes exact
two-sample crossing probabilities by lattice recursion.
"""

from __future__ import annotations

import argparse
import math


def logaddexp(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return -math.inf
    offset = max(finite)
    return offset + math.log(sum(math.exp(value - offset) for value in finite))


def log_binom(n: int, k: int) -> float:
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def log_smirnov_one_sided_tail(n: int, epsilon: float) -> float:
    """Return log P(D_n^- > epsilon) for 0 < epsilon < 1."""
    if n <= 0:
        raise ValueError("n must be positive")
    if not 0.0 < epsilon < 1.0:
        raise ValueError("epsilon must be in (0, 1)")

    terms = [n * math.log1p(-epsilon)]
    j_max = math.floor(n * (1.0 - epsilon))
    for j in range(1, j_max + 1):
        right_gap = 1.0 - epsilon - j / n
        if right_gap <= 0.0:
            terms.append(-math.inf)
            continue
        term = (
            math.log(epsilon)
            + log_binom(n, j)
            + (j - 1) * math.log(epsilon + j / n)
            + (n - j) * math.log(right_gap)
        )
        terms.append(term)
    return logaddexp(terms)


def smirnov_one_sided_tail(n: int, epsilon: float) -> float:
    return math.exp(log_smirnov_one_sided_tail(n, epsilon))


def dkwm_bound(n: int, epsilon: float) -> float:
    return min(1.0, 2.0 * math.exp(-2.0 * n * epsilon * epsilon))


def two_sample_one_sided_tail(n: int, m: int, delta: float) -> float:
    """Return P(D_{n,m}^+ >= delta) by exact lattice recursion."""
    if n <= 0 or m <= 0:
        raise ValueError("n and m must be positive")
    if delta <= 0.0:
        return 1.0

    counts = [[0] * (m + 1) for _ in range(n + 1)]
    counts[0][0] = 1
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            if i / n - j / m >= delta:
                counts[i][j] = 0
                continue
            left = counts[i - 1][j] if i > 0 else 0
            below = counts[i][j - 1] if j > 0 else 0
            counts[i][j] = left + below
    total = math.comb(n + m, n)
    return 1.0 - counts[n][m] / total


def two_sample_two_sided_tail(n: int, m: int, delta: float) -> float:
    """Return P(D_{n,m} >= delta) by exact two-wall lattice recursion."""
    if n <= 0 or m <= 0:
        raise ValueError("n and m must be positive")
    if delta <= 0.0:
        return 1.0

    counts = [[0] * (m + 1) for _ in range(n + 1)]
    counts[0][0] = 1
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            if abs(i / n - j / m) >= delta:
                counts[i][j] = 0
                continue
            left = counts[i - 1][j] if i > 0 else 0
            below = counts[i][j - 1] if j > 0 else 0
            counts[i][j] = left + below
    total = math.comb(n + m, n)
    return 1.0 - counts[n][m] / total


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--m", type=int, default=100)
    parser.add_argument(
        "--epsilon",
        type=float,
        nargs="+",
        default=[0.08, 0.10, 0.12, 0.15],
    )
    parser.add_argument("--two-sample-delta", type=float, default=0.12)
    args = parser.parse_args()

    print(f"n={args.n}")
    print("epsilon,one_sided_exact,two_sided_dkwm_bound")
    for epsilon in args.epsilon:
        exact = smirnov_one_sided_tail(args.n, epsilon)
        bound = dkwm_bound(args.n, epsilon)
        print(f"{epsilon:.6g},{exact:.12g},{bound:.12g}")
    two_sample = two_sample_one_sided_tail(args.n, args.m, args.two_sample_delta)
    two_sided = two_sample_two_sided_tail(args.n, args.m, args.two_sample_delta)
    print(
        "two_sample_one_sided:"
        f" n={args.n},m={args.m},delta={args.two_sample_delta},"
        f"tail={two_sample:.12g}"
    )
    print(
        "two_sample_two_sided:"
        f" n={args.n},m={args.m},delta={args.two_sample_delta},"
        f"tail={two_sided:.12g}"
    )


if __name__ == "__main__":
    main()
