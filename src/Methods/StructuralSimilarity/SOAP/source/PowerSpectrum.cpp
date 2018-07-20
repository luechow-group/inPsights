//
// Created by Michael Heuer on 15.05.18.
//

#include "PowerSpectrum.h"
#include "NeighborhoodExpansion.h"

//function to calculate p_ab(X_i)
Eigen::VectorXcd PowerSpectrum::partialPowerSpectrum(const NeighborhoodExpansion& n1a,
                                                     const NeighborhoodExpansion& n1b) {

    const auto & nmax = ExpansionSettings::Radial::nmax;
    const auto & lmax = ExpansionSettings::Angular::lmax;
    unsigned angularEntityLength = 2 * lmax + 1;
    unsigned entityLength = nmax * nmax * angularEntityLength;
    Eigen::VectorXcd expansionCoefficients = Eigen::VectorXcd::Zero(entityLength);

    for (unsigned n1 = 1; n1 <= nmax; ++n1) {
        unsigned n1BlockStartIdx = (n1 - 1) * nmax * angularEntityLength;

        for (unsigned n2 = 1; n2 <= nmax; ++n2) {
            unsigned n1n2BlockStartIdx = n1BlockStartIdx + (n2 - 1) * angularEntityLength;

            for (unsigned l = 0; l <= lmax; ++l) {
                //TODO NORM HERE?
                //expansionCoefficients[n1n2BlockStartIdx + l] = std::norm(powerSpectrumCoefficient(n1a, n1b,n1, n2, l));
                expansionCoefficients[n1n2BlockStartIdx + l] = powerSpectrumCoefficient(n1a, n1b,n1, n2, l);
            }
        }
    }
    return expansionCoefficients;
}

std::complex<double> PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& generic,
                                                             unsigned n1, unsigned n2, unsigned l) {
    return powerSpectrumCoefficient(generic,generic,n1,n2,l); // TODO CAN WE SAVE EFFORT HERE?
}

#include <iostream>
//double
std::complex<double> PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& speciesA,
                                                             const NeighborhoodExpansion& speciesB,
                                                             unsigned n1, unsigned n2, unsigned l ) {
    ExpansionSettings::Radial::checkBounds(n1);
    ExpansionSettings::Radial::checkBounds(n2);
    ExpansionSettings::Angular::checkBounds(l);

    std::complex<double> sum2 = (speciesA.getCoefficients_nl(n1, l).array().conjugate()
                  * speciesB.getCoefficients_nl(n2, l).array()).sum();
    //std::cout << speciesA.getCoefficients_nl(n1, l).transpose()//<< ", "<< speciesB.getCoefficients_nl(n2, l).transpose()
    //          << " ELEMENTS: ";
    //std::complex<double> sum={0,0};
    //for (int m = -int(l); m <= int(l); ++m) {
    //    //TODO sum up all particles ?? BEFORE OR AFTER MULTIPLICATION?
    //    //std::cout << speciesA.getCoefficient(n1, l, m); //<< ", " << speciesB.getCoefficient(n2, l, m);
    //    sum += std::conj(speciesA.getCoefficient(n1, l, m)) * speciesB.getCoefficient(n2, l, m); // TODO CAREFUL: the conjugation order is flipped
    //}
    //std::cout << std::endl;
    //std::cout << sum << ",\t" << sum2 << std::endl;
    return M_PI * sqrt(8./(2.*l+1)) * sum2; // in the soapxx code, this prefactor can be replaced by 1 (using a bool)
    //return M_PI * sqrt(8./(2.*l+1)) * std::norm(sum);
}

