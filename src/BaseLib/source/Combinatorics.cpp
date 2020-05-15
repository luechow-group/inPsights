/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
