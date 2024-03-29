// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "LocalSimilarity.h"
#include "ParticleKit.h"
#include "PowerSpectrum.h"
#include "NeighborhoodExpander.h"

namespace SOAP {
    namespace LocalSimilarity {
        double kernel(const Environment &e1, const Environment &e2, double zeta) {
            assert(zeta > 0 && "Zeta must be positive.");
            NeighborhoodExpander expander;
            auto expansion1 = expander.computeParticularExpansions(e1);
            auto expansion2 = expander.computeParticularExpansions(e2);

            return kernel(expansion1, expansion2, zeta);
        }

        double unnormalizedKernel(const Environment &e1, const Environment &e2) {
            NeighborhoodExpander expander;
            auto expansion1 = expander.computeParticularExpansions(e1);
            auto expansion2 = expander.computeParticularExpansions(e2);

            return unnormalizedKernel(expansion1, expansion2);
        }

        double unnormalizedSelfKernel(const Environment &e) {
            NeighborhoodExpander expander;
            auto expansion = expander.computeParticularExpansions(e);

            return unnormalizedSelfKernel(expansion);
        }

        double kernel(
                const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                const TypeSpecificNeighborhoodsAtOneCenter &expansions2, double zeta) {
            assert(zeta > 0 && "Zeta must be positive.");

            auto selfSimilarity1 = unnormalizedSelfKernel(expansions1);
            auto selfSimilarity2 = unnormalizedSelfKernel(expansions2);

            auto eps = std::numeric_limits<double>::epsilon();
            if (selfSimilarity1 <= eps || selfSimilarity2 <= eps) {
                if (selfSimilarity1 <= eps && selfSimilarity2 <= eps) {
                    spdlog::debug("LocalSelfSimilarity: Warning: "
                                  "The analyzed structure contains two isolated environments.");
                    return 1.0;
                } else {
                    spdlog::debug("LocalSelfSimilarity: Warning: "
                                  "The analyzed structure contains one isolated environments.");
                    return 0.0;
                }
            }
            auto similarityValue = unnormalizedKernel(expansions1, expansions2);
            return pow(similarityValue / sqrt(selfSimilarity1 * selfSimilarity2), zeta);
        }

        double unnormalizedKernel(
                const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {

            double similarityValue = 0;
            auto mode = General::settings.mode();
            switch (mode) {
                case General::Mode::typeAgnostic: {
                    similarityValue = internal::typeAgnostic(expansions1, expansions2);
                    break;
                }
                case General::Mode::chemical: {
                    similarityValue = internal::chemical(expansions1, expansions2);
                    break;
                }
                case General::Mode::alchemical: {
                    similarityValue = internal::alchemical(expansions1, expansions2);
                    break;
                }
                default:
                    throw std::exception();
            }
            assert(similarityValue >= 0 && "The similarity cannot be negative. "
                                           "(It might be zero if one of the centers is completely isolated)");
            return similarityValue;
        }

        double unnormalizedSelfKernel(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            double similarityValue = 0;
            switch (General::settings.mode()) {
                case General::Mode::typeAgnostic: {
                    similarityValue = internal::typeAgnostic(expansions);
                    break;
                }
                case General::Mode::chemical: {
                    similarityValue = internal::chemical(expansions);
                    break;
                }
                case General::Mode::alchemical: {
                    similarityValue = internal::alchemical(expansions, expansions);
                    // a dedicated self-similarity method for the alchemical expansion does not increase efficiency here
                    break;
                }
                default:
                    throw std::exception();
            }
            assert(similarityValue >= 0 && "The similarity cannot be negative. "
                                           "(It might be zero if the center is completely isolated)");
            return similarityValue;
        }

        double kernelDistance(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                              const TypeSpecificNeighborhoodsAtOneCenter &expansions2, double zeta) {

            return sqrt(2.0 - 2.0 * kernel(expansions1, expansions2, zeta));
        }

        double kernelDistance(const Environment &e1, const Environment &e2, double zeta) {
            assert(zeta > 0 && "Zeta must be positive.");
            NeighborhoodExpander expander;
            auto expansion1 = expander.computeParticularExpansions(e1);
            auto expansion2 = expander.computeParticularExpansions(e2);

            return kernelDistance(expansion1, expansion2, zeta);
        }



        double internal::typeAgnostic(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            int noneType = 0;
            const auto &exp = expansions.find(noneType)->second;

            auto ps = PowerSpectrum::partialPowerSpectrum(exp, exp).normalized(); // TODO CHECK THIS!
            return std::norm(ps.dot(ps)); // TODO NORM HERE OR LATER?
        }

        double internal::typeAgnostic(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                                      const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {

            int noneType = 0;
            const auto &exp1 = expansions1.find(noneType)->second;
            const auto &exp2 = expansions2.find(noneType)->second;

            auto ps1 = PowerSpectrum::partialPowerSpectrum(exp1, exp1).normalized();
            auto ps2 = PowerSpectrum::partialPowerSpectrum(exp2, exp2).normalized();
            return std::norm(ps1.dot(ps2));
        }

        double internal::chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions) {
            std::complex<double> sum = {0, 0};

            for (auto &alpha : ParticleKit::kit) {
                const auto &alphaExpansion = expansions.find(alpha.first)->second;

                for (auto &beta : ParticleKit::kit) {
                    const auto &betaExpansion = expansions.find(beta.first)->second;

                    auto ps = PowerSpectrum::partialPowerSpectrum(alphaExpansion, betaExpansion);
                    sum += ps.dot(ps);// TODO NORM HERE?
                }
            }
            return std::norm(sum);// TODO NORM HERE?
        }

        double internal::chemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                                  const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {
            std::complex<double> sum = 0;
            for (auto &alpha : ParticleKit::kit) {
                const auto &alphaExpansion1 = expansions1.find(alpha.first)->second;
                const auto &alphaExpansion2 = expansions2.find(alpha.first)->second;

                for (auto &beta : ParticleKit::kit) {
                    const auto &betaExpansion1 = expansions1.find(beta.first)->second;
                    const auto &betaExpansion2 = expansions2.find(beta.first)->second;

                    auto ps1 = PowerSpectrum::partialPowerSpectrum(alphaExpansion1, betaExpansion1);
                    auto ps2 = PowerSpectrum::partialPowerSpectrum(alphaExpansion2, betaExpansion2);
                    sum += ps1.dot(ps2);
                }
            }
            return std::norm(sum);
        }

        double internal::alchemicalTypeSimilarity(int typeA, int typeB) {
            if (typeA == typeB) {
                const auto itSame = General::settings.pairSimilarities.find({typeA, typeA});
                return (itSame != General::settings.pairSimilarities.end())? itSame->second : 1.0;
            } else {
                const auto itDiff = typeA < typeB ?
                                    General::settings.pairSimilarities.find({typeA, typeB}) :// (typeA < typeB)
                                    General::settings.pairSimilarities.find({typeB, typeA}); // (typeA > typeB)
                return (itDiff != General::settings.pairSimilarities.end()) ? itDiff->second : 0.0;
            }
        }

        double internal::alchemical(const TypeSpecificNeighborhoodsAtOneCenter &expansions1,
                                    const TypeSpecificNeighborhoodsAtOneCenter &expansions2) {
            std::complex<double> sum = 0;
            for (auto &alpha : ParticleKit::kit) {
                for (auto &alphaPrimed : ParticleKit::kit) {

                    double kappa_aap = internal::alchemicalTypeSimilarity(alpha.first, alphaPrimed.first);
                    if (kappa_aap > 0.0) {
                        const auto &alphaExpansion1 = expansions1.find(alpha.first)->second;
                        const auto &alphaPrimedExpansion2 = expansions2.find(alphaPrimed.first)->second;

                        for (auto &beta : ParticleKit::kit) {
                            for (auto &betaPrimed : ParticleKit::kit) {

                                double kappa_bbp = internal::alchemicalTypeSimilarity(beta.first, betaPrimed.first);
                                if (kappa_bbp > 0.0) {
                                    const auto &betaExpansion1 = expansions1.find(beta.first)->second;
                                    const auto &betaPrimedExpansion2 = expansions2.find(betaPrimed.first)->second;

                                    auto ps1 = PowerSpectrum::partialPowerSpectrum(alphaExpansion1, betaExpansion1);
                                    auto ps2 = PowerSpectrum::partialPowerSpectrum(alphaPrimedExpansion2,
                                                                                   betaPrimedExpansion2);

                                    sum += ps1.dot(ps2) * kappa_aap * kappa_bbp;
                                }
                            }
                        }
                    }
                }
            }
            return std::norm(sum);
        }
    }
}
