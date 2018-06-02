//
// Created by Michael Heuer on 15.05.18.
//

#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"

double LocalSimilarity::localSimilarity(const Environment& e1, const Environment& e2, double zeta) {
    assert(zeta > 0 && "Zeta must be positive.");
    NeighborhoodExpander expander;
    auto exp1 = expander.computeParticularExpansions(e1);
    auto exp2 = expander.computeParticularExpansions(e2);

    return pow(localSimilarity(exp1,exp2),zeta);
}

double LocalSimilarity::unnormalizedLocalSimialrity(const Environment& e1,
                                                    const Environment& e2) {
    NeighborhoodExpander expander;
    auto exp1 = expander.computeParticularExpansions(e1);
    auto exp2 = expander.computeParticularExpansions(e2);

    return unnormalizedLocalSimilarity(exp1, exp2);
}

double LocalSimilarity::localSimilarity(
        const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
        const TypeSpecificNeighborhoodsAtOneCenter& expansions2, double zeta) {
    assert(zeta > 0 && "Zeta must be positive.");

    auto similarityValue = unnormalizedLocalSimilarity(expansions1, expansions2)
                           / sqrt(unnormalizedLocalSelfSimilarity(expansions1)
                                 *unnormalizedLocalSelfSimilarity(expansions2));
    /// sqrt(unnormalizedLocalSimilarity(expansions1,expansions1)
    //      *unnormalizedLocalSimilarity(expansions2,expansions2));


    return pow(similarityValue,zeta);
}

double LocalSimilarity::unnormalizedLocalSimilarity(
        const TypeSpecificNeighborhoodsAtOneCenter& expansions1,
        const TypeSpecificNeighborhoodsAtOneCenter& expansions2) { //is type information of the neighbors needed here?
    //TODO MODIFY

    double similarityValue = 0;
    switch (ExpansionSettings::mode) {
        case ExpansionSettings::Mode::Generic: {

            int noneType = 0;
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
        case ExpansionSettings::Mode::TypeSpecific: {

            double sumAB = 0;

            for (auto & typeA : ParticleKit::kit) {
                const auto &e1a = expansions1.find(typeA.first)->second;
                const auto &e2a = expansions2.find(typeA.first)->second;
                //std::cout << e1a << std::endl<< std::endl;
                //std::cout << e2a << std::endl<< std::endl;

                for (auto & typeB : ParticleKit::kit) {
                    // Alchemical similarity can be implemented here
                    if (typeA.first == typeB.first) { //TODO kronecker delta or
                        const auto &e1b = expansions1.find(typeB.first)->second;
                        const auto &e2b = expansions2.find(typeB.first)->second;
                        //std::cout << e1b << std::endl<< std::endl;
                        //std::cout << e2b << std::endl<< std::endl;

                        auto ps1 = PowerSpectrum::partialPowerSpectrum(e1a, e1b); // kroneckerdelta => (e1a, e1a)
                        auto ps2 = PowerSpectrum::partialPowerSpectrum(e2a, e2b); // kroneckerdelta => (e2a, e2a)
                        //std::cout << ps1.transpose() << std::endl;
                        //std::cout << ps2.transpose() << std::endl;

                        sumAB += ps1.dot(ps2);
                    }
                }
                /*auto ps1 = PowerSpectrum::partialPowerSpectrum(e1a, e1a); // kroneckerdelta => (e1a, e1a)
                auto ps2 = PowerSpectrum::partialPowerSpectrum(e2a, e2a); // kroneckerdelta => (e2a, e2a)
                sumAB += ps1.dot(ps2);*/
            }
            similarityValue = sumAB;
        }
    }
    return similarityValue;
}

double LocalSimilarity::unnormalizedLocalSelfSimilarity(const TypeSpecificNeighborhoodsAtOneCenter& expansions) {
    double similarityValue = 0;
    switch (ExpansionSettings::mode) {
        case ExpansionSettings::Mode::Generic: {
            int noneType = 0;
            const auto &expansion = expansions.find(noneType)->second;
            auto partialPowerSpectrum = PowerSpectrum::partialPowerSpectrum(expansion, expansion).normalized();

            similarityValue = partialPowerSpectrum.dot(partialPowerSpectrum);
            break;
        }
        case ExpansionSettings::Mode::TypeSpecific: {

            double sumAB = 0;

            for (auto & typeA : ParticleKit::kit) {
                const auto &expansionA = expansions.find(typeA.first)->second;

                for (auto & typeB : ParticleKit::kit) {
                    // Alchemical similarity can be implemented here
                    if (typeA.first == typeB.first) { //TODO kronecker delta
                        const auto &expansionB = expansions.find(typeB.first)->second;

                        auto partialPowerSpectrum = PowerSpectrum::partialPowerSpectrum(expansionA, expansionB); // kroneckerdelta => (e1a, e1a)
                        sumAB += partialPowerSpectrum.dot(partialPowerSpectrum);
                    }
                }
            }
            similarityValue = sumAB;
        }
    }
    return similarityValue;
}

/*static std::complex<double> distance(
        const PowerSpectrum &a,
        const PowerSpectrum &b) {
    return sqrt(2 + 2 * a.asEigenVector().normalized().dot(b.asEigenVector().normalized()));
}*/
