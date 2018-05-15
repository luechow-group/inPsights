//
// Created by Michael Heuer on 09.05.18.
//

#ifndef AMOLQCPP_LOCALSIMILARITY_H
#define AMOLQCPP_LOCALSIMILARITY_H

#include <complex>

class Environment;

namespace LocalSimilarity {

    double unnormalizedLocalSimialrity(const Environment& e1, const Environment& e2);
    double localSimilarity(const Environment& e1, const Environment& e2);

    /*static std::complex<double> distance(
            const PowerSpectrum &a,
            const PowerSpectrum &b) {
        return sqrt(2 + 2 * a.asEigenVector().normalized().dot(b.asEigenVector().normalized()));
    }*/

};

#endif //AMOLQCPP_LOCALSIMILARITY_H
