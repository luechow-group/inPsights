//
// Created by Michael Heuer on 15.05.18.
//

#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"


namespace LocalSimilarity {
    double localSimilarity(const Environment &e1, const Environment &e2, double zeta) {
        assert(zeta > 0 && "Zeta must be positive.");
        NeighborhoodExpander expander;
        auto exp1 = expander.computeParticularExpansions(e1);
        auto exp2 = expander.computeParticularExpansions(e2);

        return pow(localSimilarity(exp1, exp2), zeta);
    }

    double unnormalizedLocalSimialrity(const Environment &e1,
                                                        const Environment &e2) {
        NeighborhoodExpander expander;
        auto exp1 = expander.computeParticularExpansions(e1);
        auto exp2 = expander.computeParticularExpansions(e2);

        return unnormalizedLocalSimilarity(exp1, exp2);
    }

    double localSimilarity(
            const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
            const TypeSpecificNeighborhoodsAtOneCenter &expansions2, double zeta) {
        assert(zeta > 0 && "Zeta must be positive.");

        auto similarityValue = unnormalizedLocalSimilarity(expansions1, expansions2)
                               / sqrt(unnormalizedLocalSelfSimilarity(expansions1)
                                      *unnormalizedLocalSelfSimilarity(expansions2));
        /// sqrt(unnormalizedLocalSimilarity(expansions1,expansions1)
        //      *unnormalizedLocalSimilarity(expansions2,expansions2));

        return pow(similarityValue, zeta);
    }

    double unnormalizedLocalSimilarity(
            const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
            const TypeSpecificNeighborhoodsAtOneCenter &expansions2) { //is type information of the neighbors needed here?

        double similarityValue = 0;
        switch (ExpansionSettings::mode) {
            case ExpansionSettings::Mode::Generic: {
                similarityValue = generic(expansions1, expansions2);
                break;
            }
            case ExpansionSettings::Mode::Chemical: {
                similarityValue = typeSpecific(expansions1, expansions2);
                break;
            }
            case ExpansionSettings::Mode::Alchemical: {
                similarityValue = alchemical(expansions1,expansions2);
                break;
            }
        }
        return similarityValue;
    }

    // use factory method and eliminate switch
    double unnormalizedLocalSelfSimilarity(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
        double similarityValue = 0;
        switch (ExpansionSettings::mode) {
            case ExpansionSettings::Mode::Generic: {
                similarityValue = generic(expansions);
                break;
            }
            case ExpansionSettings::Mode::Chemical: {
                similarityValue = typeSpecific(expansions);
                break;
            }
            case ExpansionSettings::Mode::Alchemical: {
                similarityValue = alchemical(expansions,expansions); //TODO: self similarity
                break;
            }
        }
        return similarityValue;
    }

    namespace {
        double generic(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            int noneType = 0;
            const auto &exp = expansions.find(noneType)->second;
            auto ps = PowerSpectrum::partialPowerSpectrum(exp, exp).normalized();

            return ps.normalized().dot(ps);
        }

        double generic(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                       const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {

            int noneType = 0;
            const auto &exp1 = expansions1.find(noneType)->second;
            const auto &exp2 = expansions2.find(noneType)->second;

            return PowerSpectrum::partialPowerSpectrum(exp1, exp1).normalized().dot(
                    PowerSpectrum::partialPowerSpectrum(exp2, exp2).normalized());
        }

        double typeSpecific(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            double sum = 0;

            for (auto &alpha : ParticleKit::kit) {
                const auto &alphaExpansion = expansions.find(alpha.first)->second;

                for (auto &beta : ParticleKit::kit) {
                    // Alchemical similarity can be implemented here
                    //if (alpha.first == beta.first) { //TODO kronecker delta?
                        const auto &betaExpansion = expansions.find(beta.first)->second;

                        // TODO NORMALIZATION CORRECT? PAPER USES UNNORMALIZED PARTIAL POWER SPECTRA
                        auto partialPowerSpectrum = PowerSpectrum::partialPowerSpectrum(alphaExpansion, betaExpansion).normalized();
                        sum += partialPowerSpectrum.dot(partialPowerSpectrum);
                    //}
                }
            }
            return sum;
        }

        double typeSpecific(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                            const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {
            double sum = 0;
            for (auto &alpha : ParticleKit::kit) {
                const auto &alphaExpansion1 = expansions1.find(alpha.first)->second;
                const auto &alphaExpansion2 = expansions2.find(alpha.first)->second;

                for (auto &beta : ParticleKit::kit) {
                    //if (alpha.first == beta.first) { //TODO kronecker delta?
                        const auto &betaExpansion1 = expansions1.find(beta.first)->second;
                        const auto &betaExpansion2 = expansions2.find(beta.first)->second;

                        // TODO NORMALIZATION CORRECT? PAPER USES UNNORMALIZED PARTIAL POWER SPECTRA
                        auto ps1 = PowerSpectrum::partialPowerSpectrum(alphaExpansion1, betaExpansion1).normalized();
                        auto ps2 = PowerSpectrum::partialPowerSpectrum(alphaExpansion2, betaExpansion2).normalized();

                        sum += ps1.dot(ps2);
                    //}
                }
            }
            return sum;
        }

        double kroneckerDelta(int typeA, int typeB) {
            std::map<std::pair<int, int>, double>::const_iterator it;

            if (typeA == typeB)
                return 1.0;
            else if (typeA < typeB)
                it = ExpansionSettings::Alchemical::pairSimilarities.find({typeA, typeB});
            else
                it = ExpansionSettings::Alchemical::pairSimilarities.find({typeB, typeA});


            if (it != ExpansionSettings::Alchemical::pairSimilarities.end())
                return (*it).second;
            else
                return 0.0;
        }

        /* TODO
        double alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            double sum = 0;
            for (auto &alpha : ParticleKit::kit) {
                for (auto &alphaPrimed : ParticleKit::kit) {

                    double k_aap = kroneckerDelta(alpha.first, alphaPrimed.first);
                    if (k_aap > 0.0) {
                        const auto &alphaExpansion = expansions.find(alpha.first)->second;

                        for (auto &beta : ParticleKit::kit) {
                            for (auto &betaPrimed : ParticleKit::kit) {

                                double k_bbp = kroneckerDelta(beta.first, betaPrimed.first);
                                if (k_bbp> 0.0) {
                                    const auto &betaExpansion = expansions.find(beta.first)->second;

                                    auto ps = PowerSpectrum::partialPowerSpectrum(alphaExpansion,
                                                                                   betaExpansion);

                                    sum += ps.dot(ps) * k_aap*k_bbp;
                                }
                            }
                        }
                    }
                }
            }
            return sum;
        }*/

        double alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                          const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {
            double sum = 0;
            for (auto &alpha : ParticleKit::kit) {
                for (auto &alphaPrimed : ParticleKit::kit) {

                    double k_aap = kroneckerDelta(alpha.first, alphaPrimed.first);
                    if (k_aap > 0.0) {
                        const auto &alphaExpansion1 = expansions1.find(alpha.first)->second;
                        const auto &alphaPrimedExpansion2 = expansions2.find(alphaPrimed.first)->second;
                        //TODO WHAT IF NOT FOUND? it == end()

                        for (auto &beta : ParticleKit::kit) {
                            for (auto &betaPrimed : ParticleKit::kit) {

                                double k_bbp = kroneckerDelta(beta.first, betaPrimed.first);
                                if (k_bbp > 0.0) {
                                    const auto &betaExpansion1 = expansions1.find(beta.first)->second;
                                    const auto &betaPrimedExpansion2 = expansions2.find(betaPrimed.first)->second;

                                    auto ps1 = PowerSpectrum::partialPowerSpectrum(alphaExpansion1,
                                                                                   betaExpansion1);
                                    auto ps2 = PowerSpectrum::partialPowerSpectrum(alphaPrimedExpansion2,
                                                                                   betaPrimedExpansion2);

                                    sum += ps1.dot(ps2) * k_aap*k_bbp;
                                }
                            }
                        }
                    }
                }
            }
            return sum;
        }
    }
}
/*static std::complex<double> distance(
        const PowerSpectrum &a,
        const PowerSpectrum &b) {
    return sqrt(2 + 2 * a.asEigenVector().normalized().dot(b.asEigenVector().normalized()));
}*/
