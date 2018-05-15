//
// Created by Michael Heuer on 15.05.18.
//

#include "LocalSimilarity.h"
#include "ParticlePool.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"

#include <exception>

double LocalSimilarity::unnormalizedLocalSimialrity(const Environment& e1,
                                                    const Environment& e2) {
    NeighborhoodExpander neighborhoodExpander;
    // collect expansions
    std::map<Elements::ElementType, NeighborhoodExpansion> expansions1;
    std::map<Elements::ElementType, NeighborhoodExpansion> expansions2;

    auto numberOfElementTypes = 0;

    switch (ExpansionSettings::mode){
        case ExpansionMode::Generic:{

            auto noneType = Elements::ElementType::none;
            auto exp1 = neighborhoodExpander.expandEnvironment(e1, noneType);
            auto exp2 = neighborhoodExpander.expandEnvironment(e2, noneType);

            std::cout << "exp1\n" << exp1 << std::endl;//TODO Delete
            std::cout << "exp2\n" << exp2 << std::endl;

            auto ps1 = PowerSpectrum::partialPowerSpectrum(exp1, exp1).normalized();
            auto ps2 = PowerSpectrum::partialPowerSpectrum(exp2, exp2).normalized();

            std::cout << "ps1\n" << ps1.transpose() << std::endl;
            std::cout << "ps2\n" << ps2.transpose()  << std::endl;

            return ps1.dot(ps2);
        }

        case ExpansionMode::TypeSpecific:{

            numberOfElementTypes = unsigned(ParticlePool::atomKit().size());

            for (unsigned t = 0; t < numberOfElementTypes; ++t) {
                auto elementType = ParticlePool::atomKit()[t].first;
                expansions1.insert({elementType, neighborhoodExpander.expandEnvironment(e1, elementType)});
                expansions2.insert({elementType, neighborhoodExpander.expandEnvironment(e2, elementType)});
            };

            double sumAB = 0;

            for (unsigned a = 0; a < numberOfElementTypes; ++a) {

                auto typeA =ParticlePool::atomKit()[a].first;

                auto exp1a = expansions1.find(typeA)->second;
                auto exp2a = expansions2.find(typeA)->second;

                std::cout << exp1a << std::endl;//TODO Delete
                std::cout << exp2a << std::endl;

                for (unsigned b = 0; b < numberOfElementTypes; ++b) {
                    //p is the powerspectrum class
                    auto typeB = ParticlePool::atomKit()[b].first;

                    if(typeA == typeB) { // kronecker delta
                        auto exp1b = expansions1.find(typeB)->second;
                        auto exp2b = expansions2.find(typeB)->second;

                        std::cout << exp1b << std::endl;//TODO Delete
                        std::cout << exp2b << std::endl;

                        auto ps1 = PowerSpectrum::partialPowerSpectrum(exp1a, exp1b);
                        auto ps2 = PowerSpectrum::partialPowerSpectrum(exp2a, exp2b);

                        if( ExpansionSettings::mode == ExpansionMode::Generic){
                            ps1.normalize();
                            ps2.normalize();
                        }

                        sumAB += ps1.dot(ps2);
                    }
                }
            }

            return sumAB;
        }
    }
}


/*static std::complex<double> distance(
        const PowerSpectrum &a,
        const PowerSpectrum &b) {
    return sqrt(2 + 2 * a.asEigenVector().normalized().dot(b.asEigenVector().normalized()));
}*/

//atom specific
double LocalSimilarity::localSimilarity(const Environment& e1, const Environment& e2) {
    return unnormalizedLocalSimialrity(e1, e2)
           / sqrt(unnormalizedLocalSimialrity(e1,e1) * unnormalizedLocalSimialrity(e2,e2));
}