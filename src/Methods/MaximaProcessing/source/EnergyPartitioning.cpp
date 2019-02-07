//
// Created by Michael Heuer on 2019-02-03.
//

#include <EnergyPartitioning.h>
#include <algorithm>
#include <Motif.h>

namespace EnergyPartitioning {
    namespace MotifBased{
        MotifEnergies calculateInterationEnergies(const Motifs &motifs,
                                                  const EnergyStatistics::ElectronicEnergy &electronicEnergy,
                                                  const Eigen::MatrixXd &Vnn) {
            MotifEnergies motifEnergies{};

            for(auto motif = motifs.motifVector_.begin(); motif != motifs.motifVector_.end(); ++motif) {

                // self-interaction
                auto intraEnergy = calculateSelfInteractionEnergy(*motif, electronicEnergy, Vnn);
                motifEnergies.addPair({*motif, *motif}, intraEnergy);

                // interaction with other motifs
                std::map<Eigen::Index, double> interEnergies;
                for (auto otherMotif = motif; otherMotif != motifs.motifVector_.end(); ++otherMotif) {
                    auto interEnergy = caclulateInteractionEnergy(*motif, *otherMotif, electronicEnergy, Vnn);
                    motifEnergies.addPair({*motif, *otherMotif}, interEnergy);
                }
            }
            return motifEnergies;
        }

        double caclulateInteractionEnergy(const Motif &motif, const Motif &otherMotif,
                                          const EnergyStatistics::ElectronicEnergy &electronicEnergy,
                                          const Eigen::MatrixXd &Vnn) {
            double inter = 0;

            for (auto i : motif.electronIndices())
                for (auto j : otherMotif.electronIndices())
                    inter += electronicEnergy.Vee().mean()(i, j); // Vee


            for (auto i : motif.electronIndices())
                for (auto l : otherMotif.atomIndices()) {

                    inter += electronicEnergy.Ven().mean()(i, l); // Ven ->
                }

            for (auto j : otherMotif.electronIndices())
                for (auto k : motif.atomIndices())
                    inter += electronicEnergy.Ven().mean()(j, k); // Ven <-


            for (auto k : motif.atomIndices())
                for (auto l : otherMotif.atomIndices())
                    inter +=  Vnn(k, l); // Vnn

            return inter;
        }

        double calculateSelfInteractionEnergy(
                const Motif &motif,
                const EnergyStatistics::ElectronicEnergy &electronicEnergy,
                const Eigen::MatrixXd &Vnn
                ) {

            double intra = 0;

            for (auto i = motif.electronIndices().begin(); i != motif.electronIndices().end(); ++i) {
                intra += electronicEnergy.Te().mean()(*i); // Te

                for (auto k = motif.atomIndices().begin(); k != motif.atomIndices().end(); ++k)
                    intra += electronicEnergy.Ven().mean()(*i, *k); // Ven

            }

            for (auto i = motif.electronIndices().begin(); i != std::prev(motif.electronIndices().end()); ++i) {
                for (auto j = next(i); j != motif.electronIndices().end(); ++j)
                    intra += electronicEnergy.Vee().mean()(*i, *j); // Vee
            }

            // atoms in motif ( zero contribution if zero or one atom is present)
            for (auto k = motif.atomIndices().begin(); k != std::prev(motif.atomIndices().end()); ++k)
                for (auto l = std::next(motif.atomIndices().begin()); l != motif.atomIndices().end(); ++l)
                    intra += Vnn(*k, *l); // Vnn

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
