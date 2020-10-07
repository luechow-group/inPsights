// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include <Combinatorics.h>

std::size_t Combinatorics::binomial(std::size_t n, std::size_t k) {
    assert( n <= MAX_BINOMIAL
    && "The maximal binomial fitting into an unsigned long int with all of its possible k values is n=34.");
    if (k > n) {
        return 0;
    }
    std::size_t r = 1;
    for (std::size_t d = 1; d <= k; ++d) {
        r *= n--;
        r /= d;
    }
    return r;
}

std::size_t Combinatorics::factorial(std::size_t n) {
    assert(n <=  MAX_FACTORIAL
    && "The maximal factorial fitting into an unsigned long int is n=12.");
    return (n == 0 || n == 1) ? 1 : factorial(n - 1) * n;
}
