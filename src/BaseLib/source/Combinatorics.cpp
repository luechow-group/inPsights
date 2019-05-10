//
// Created by heuer on 08.05.19.
//

#include <Combinatorics.h>

std::size_t Combinatorics::binomial(std::size_t n, std::size_t k) {
    assert( n <= 34
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
    assert(n <=  12
    && "The maximal factorial fitting into an unsigned long int is n=12.");
    return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
