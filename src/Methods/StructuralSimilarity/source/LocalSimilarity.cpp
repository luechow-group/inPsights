//
// Created by Michael Heuer on 15.05.18.
//

#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"

#include <exception>
#include <Type.h>

double LocalSimilarity::localSimilarity(const Environment& e1, const Environment& e2, unsigned zeta) {

    NeighborhoodExpander expander;
    auto exp1 = expander.computeExpansions(e1); //TODO optimization: store somewhere else
    auto exp2 = expander.computeExpansions(e2);

    auto similarityValue = unnormalizedLocalSimilarity(exp1, exp2)
            / sqrt(unnormalizedLocalSimilarity(exp1, exp1) * unnormalizedLocalSimilarity(exp2, exp2));

    std::cout << "normalized: " << similarityValue << std::endl;
    return pow(similarityValue,zeta);
}

double LocalSimilarity::unnormalizedLocalSimialrity(const Environment& e1,
                                                    const Environment& e2) {
    NeighborhoodExpander expander;
    auto exp1 = expander.computeExpansions(e1);//TODO optimization: store somewhere else
    auto exp2 = expander.computeExpansions(e2);

    auto similarityValue = unnormalizedLocalSimilarity(exp1, exp2);
    std::cout << "unnormalized:" << similarityValue << std::endl;
    return similarityValue;
}

double LocalSimilarity::unnormalizedLocalSimilarity(
        const std::map<int, NeighborhoodExpansion> &expansions1,
        const std::map<int, NeighborhoodExpansion> &expansions2) {

    //TODO MODIFY

    double similarityValue = 0;
    switch (ExpansionSettings::mode) {
        case ExpansionMode::Generic: {

            auto noneType = int(Type::None);
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
            auto numberOfTypes = ParticleKit::numberOfTypes();

            for (unsigned a = 0; a < numberOfTypes; ++a) {


                //TODO THIS IS SUPER UGLY
                int typeA;
                if(a < ParticleKit::atomKit.size()) {
                    typeA = int(ParticleKit::atomKit[a].first);//WRITE PARTICLE KIT ACCESSOR
                } else if ( ParticleKit::atomKit.size() ) {
                    typeA = int(Spins::SpinType::alpha);
                } else {//if ( ParticleKit::atomKit.size() +1 )
                    typeA = int(Spins::SpinType::beta);
                }

                const auto &e1a = expansions1.find(typeA)->second;
                const auto &e2a = expansions2.find(typeA)->second;

                for (unsigned b = 0; b < numberOfTypes; ++b) {
                    //p is the powerspectrum class
                    //TODO THIS IS SUPER UGLY
                    //TODO optimization: store differently without an if
                    int typeB;
                    if(b < ParticleKit::atomKit.size()) {
                        typeB = int(ParticleKit::atomKit[b].first);
                    } else if ( ParticleKit::atomKit.size() ) {
                        typeB = int(Spins::SpinType::alpha);
                    } else {//if ( ParticleKit::atomKit.size() +1 )
                        typeB = int(Spins::SpinType::beta);
                    }

                    // Alchemical similarity can be implemented here
                    if (typeA == typeB) { // kronecker delta
                        const auto &e1b = expansions1.find(typeB)->second;
                        const auto &e2b = expansions2.find(typeB)->second;

                        //std::cout << e1b << std::endl;//TODO Delete
                        //std::cout << e2b << std::endl;

                        auto ps1 = PowerSpectrum::partialPowerSpectrum(e1a, e1b); // kroneckerdelta => (e1a, e1a)
                        auto ps2 = PowerSpectrum::partialPowerSpectrum(e2a, e2b); // kroneckerdelta => (e2a, e2a)

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
