//
// Created by Michael Heuer on 15.05.18.
//

#include "LocalSimilarity.h"
#include "ParticlePool.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"

#include <exception>

double LocalSimilarity::localSimilarity(const Environment& e1, const Environment& e2) {

    auto exp1 = computeExpansions(e1);
    auto exp2 = computeExpansions(e2);

    auto similarityValue = unnormalizedLocalSimialrity(exp1,exp2)
            / sqrt(unnormalizedLocalSimialrity(exp1,exp1) * unnormalizedLocalSimialrity(exp2,exp2));

    std::cout << "normalized: " << similarityValue << std::endl;
    return similarityValue;
}
double LocalSimilarity::unnormalizedLocalSimialrity(const Environment& e1,
                                                    const Environment& e2) {
    auto exp1 = computeExpansions(e1);
    auto exp2 = computeExpansions(e2);

    auto similarityValue = unnormalizedLocalSimialrity(exp1,exp2);
    std::cout << "unnormalized:" << similarityValue << std::endl;
    return similarityValue;
}

std::map<Elements::ElementType, NeighborhoodExpansion>
LocalSimilarity::computeExpansions(const Environment &e) {

    NeighborhoodExpander neighborhoodExpander;
    std::map<Elements::ElementType, NeighborhoodExpansion> expansions;

    switch (ExpansionSettings::mode) {
        case ExpansionMode::Generic: {

            auto noneType = Elements::ElementType::none;
            expansions.emplace(noneType, neighborhoodExpander.expandEnvironment(e, noneType));
            break;
        }
        case ExpansionMode::TypeSpecific: {
            auto numberOfElementTypes = unsigned(ParticlePool::atomKit.size());

            for (unsigned t = 0; t < numberOfElementTypes; ++t) {
                auto elementType = ParticlePool::atomKit[t].first;
                expansions.emplace(elementType, neighborhoodExpander.expandEnvironment(e, elementType));
            };
            break;
        }
    }
    return expansions;
}

double LocalSimilarity::unnormalizedLocalSimialrity(
        const std::map<Elements::ElementType, NeighborhoodExpansion>& expansions1,
        const std::map<Elements::ElementType, NeighborhoodExpansion>& expansions2) {

    double similarityValue = 0;
    switch (ExpansionSettings::mode) {
        case ExpansionMode::Generic: {

            auto noneType = Elements::ElementType::none;
            const auto &e1 = expansions1.find(noneType)->second;
            const auto &e2 = expansions2.find(noneType)->second;

            //TODO delete
            //std::cout << e1.asEigenVector().transpose() << std::endl;
            //std::cout << e2.asEigenVector().transpose() << std::endl<< std::endl;
            auto ps1 = PowerSpectrum::partialPowerSpectrum(e1, e1).normalized();
            auto ps2 = PowerSpectrum::partialPowerSpectrum(e2, e2).normalized();

            //TODO Delete
            //std::cout << ps1.transpose() << std::endl;
            //std::cout << ps2.transpose() << std::endl<< std::endl;

            similarityValue = ps1.dot(ps2);
            break;
        }
        case ExpansionMode::TypeSpecific: {

            double sumAB = 0;
            auto numberOfElementTypes = unsigned(ParticlePool::atomKit.size());

            for (unsigned a = 0; a < numberOfElementTypes; ++a) {

                auto typeA = ParticlePool::atomKit[a].first;

                const auto &e1a = expansions1.find(typeA)->second;
                const auto &e2a = expansions2.find(typeA)->second;

                for (unsigned b = 0; b < numberOfElementTypes; ++b) {
                    //p is the powerspectrum class
                    auto typeB = ParticlePool::atomKit[b].first;

                    if (typeA == typeB) { // kronecker delta
                        const auto &e1b = expansions1.find(typeB)->second;
                        const auto &e2b = expansions2.find(typeB)->second;

                        //std::cout << e1b << std::endl;//TODO Delete
                        //std::cout << e2b << std::endl;

                        auto ps1 = PowerSpectrum::partialPowerSpectrum(e1a, e1b);
                        auto ps2 = PowerSpectrum::partialPowerSpectrum(e2a, e2b);

                        sumAB += ps1.dot(ps2);
                    }
                }
            }
            similarityValue = sumAB;
        }
    }
    std::cout << "\tunnormalized:" << similarityValue << std::endl;
    return similarityValue;
}

/*static std::complex<double> distance(
        const PowerSpectrum &a,
        const PowerSpectrum &b) {
    return sqrt(2 + 2 * a.asEigenVector().normalized().dot(b.asEigenVector().normalized()));
}*/

//atom specific

// TODO STORE
