//
// Created by Michael Heuer on 15.05.18.
//

#include "PowerSpectrum.h"
#include "NeighborhoodExpansion.h"

//function to calculate p_ab(X_i)
Eigen::VectorXd PowerSpectrum::partialPowerSpectrum(const NeighborhoodExpansion& n1a,
                                                    const NeighborhoodExpansion& n1b) {

    //TODO? check if subset with particlepool

    const auto & nmax = ExpansionSettings::Radial::nmax;
    const auto & lmax = ExpansionSettings::Angular::lmax;
    unsigned angularEntityLength = 2 * lmax + 1;
    unsigned entityLength = nmax * nmax * angularEntityLength;
    Eigen::VectorXd expansionCoefficients = Eigen::VectorXd::Zero(entityLength);

    for (unsigned n1 = 1; n1 <= nmax; ++n1) {
        unsigned n1BlockStartIdx = (n1 - 1) * nmax * angularEntityLength;

        for (unsigned n2 = 1; n2 <= nmax; ++n2) {
            unsigned n1n2BlockStartIdx = n1BlockStartIdx + (n2 - 1) * angularEntityLength;

            for (unsigned l = 0; l <= lmax; ++l) {
                //TODO NORM HERE
                //expansionCoefficients[n1n2BlockStartIdx + l] = std::norm(powerSpectrumCoefficient(n1a, n1b,n1, n2, l));
                expansionCoefficients[n1n2BlockStartIdx + l] = powerSpectrumCoefficient(n1a, n1b,n1, n2, l);
            }
        }
    }
    return expansionCoefficients;
}

double PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& generic,/*replace GenericType generic*/
                                               unsigned n1, unsigned n2, unsigned l) {
    return powerSpectrumCoefficient(generic,generic,n1,n2,l);
}

double PowerSpectrum::powerSpectrumCoefficient(const NeighborhoodExpansion& speciesA,
                                               const NeighborhoodExpansion& speciesB,
                                               unsigned n1, unsigned n2, unsigned l ) {/*, kappa_aa,kappa_b */

    //assert(speciesA.getSettings() == speciesB.getSettings());

    ExpansionSettings::Radial::checkBounds(n1);
    ExpansionSettings::Radial::checkBounds(n2);
    ExpansionSettings::Angular::checkBounds(l);

    double sum = 0;

    for (int m = -int(l); m < int(l); ++m) {
        //TODO sum up all particles ?? BEFORE OR AFTER MULTIPLICATION?
        // ACCOUNTED FOR
                auto a = speciesA.getCoefficient(n1, l, m);
                auto b = speciesA.getCoefficient(n2, l, m);
        sum += std::norm(std::conj(speciesA.getCoefficient(n1, l, m)) * speciesB.getCoefficient(n2, l, m));
        assert(sum == sum && "Sum cannot be NaN!");
    }



    return M_PI * sqrt(8./(2.*l+1)) * sum;
}

