//
// Created by Michael Heuer on 2019-02-03.
//

#include <EnergyPartitioning.h>
#include <MotifAnalysis.h>
#include <algorithm>

namespace EnergyPartitioning {

    namespace MotifBased{
        MotifEnergies calculateInterationEnergies(const MotifAnalysis::Motifs &motifs,
                                                  const EnergyStatistics::ElectronicEnergy &electronicEnergy) {
            MotifEnergies motifEnergies{};

            for(auto motif = motifs.motifVector.begin(); motif != motifs.motifVector.end(); ++motif) {
                // self-interaction
                auto intraEnergy = calculateSelfInteractionEnergy(*motif, electronicEnergy);
                motifEnergies.addPair({*motif, *motif}, intraEnergy);

                // interaction with other motifs
                std::map<Eigen::Index, double> interEnergies;
                for (auto otherMotif = motif; otherMotif != motifs.motifVector.end(); ++otherMotif) {
                    auto interEnergy = caclulateInteractionEnergy(*motif, *otherMotif, electronicEnergy);
                    motifEnergies.addPair({*motif, *otherMotif}, interEnergy);
                }
            }
            return motifEnergies;
        }

        double caclulateInteractionEnergy(const MotifAnalysis::Motif &motif, const MotifAnalysis::Motif &otherMotif,
                                          const EnergyStatistics::ElectronicEnergy &electronicEnergy) {
            double inter = 0;
            for (auto i : motif.electronIndices)
                for (auto j : otherMotif.electronIndices)
                    inter += electronicEnergy.Vee().mean()(i,j);

            return inter;
        }

        double calculateSelfInteractionEnergy(
                const MotifAnalysis::Motif &motif,
                const EnergyStatistics::ElectronicEnergy &electronicEnergy
                ) {

            double intra = 0;

            for (auto i = motif.electronIndices.begin(); i != prev(motif.electronIndices.end()); ++i) {
                intra += electronicEnergy.Te().mean()(*i)
                        + electronicEnergy.Ven().mean().row(*i).sum();

                for (auto j = next(i); j != motif.electronIndices.end(); ++j)
                    intra += electronicEnergy.Vee().mean()(*i, *j);
            }
            return intra;
        }
    }


    namespace ParticleBased {
        Eigen::VectorXd oneAtomEnergies(const Eigen::MatrixXd &Vnn) {
            auto nAtoms = Vnn.rows();

            Eigen::VectorXd En = Eigen::VectorXd::Zero(nAtoms);

            for (Eigen::Index k = 0; k < nAtoms; ++k) {
                for (Eigen::Index l = k + 1; l < nAtoms; ++l) {
                    En[k] += 0.5 * Vnn(k, l);
                    En[l] += 0.5 * Vnn(k, l);
                }
            }
            return En;
        }

        Eigen::VectorXd oneElectronEnergies(const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                            const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            auto nElectrons = Ven.rows();
            auto nAtoms = Ven.cols();

            Eigen::VectorXd Ee = Eigen::VectorXd::Zero(nElectrons);

            for (Eigen::Index i = 0; i < nElectrons; ++i) {
                Ee[i] += Te[i];

                for (Eigen::Index k = 0; k < nAtoms; ++k)
                    Ee[i] += Ven(i, k);

                for (Eigen::Index j = i + 1; j < nElectrons; ++j) {
                    Ee[i] += 0.5 * Vee(i, j);
                    Ee[j] += 0.5 * Vee(i, j);
                }
            }
            return Ee;
        }
    }
}
