//
// Created by Michael Heuer on 2019-02-03.
//

#include <EnergyPartitioning.h>
#include <algorithm>
#include <Motif.h>

namespace EnergyPartitioning {

    double calculateTotalEnergy(const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {

        double totalEnergy = 0;

        auto nElectrons = Ven.rows();
        auto nAtoms = Ven.cols();

        for (Eigen::Index i = 0; i < nElectrons; ++i) {
            totalEnergy += Te[i];

            for (Eigen::Index k = 0; k < nAtoms; ++k)
                totalEnergy += Ven(i, k);

            for (Eigen::Index j = i + 1; j < nElectrons; ++j)
                totalEnergy += Vee(i, j);
        }

        for (Eigen::Index k = 0; k < nAtoms; ++k)
            for (Eigen::Index l = k+1; l < nAtoms; ++l)
                totalEnergy += Vnn(k, l);

        return totalEnergy;
    }

    namespace MotifBased{
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> calculateInterationEnergies(const Motifs &motifs,
                                                  const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                  const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {

            Eigen::VectorXd intraEnergies = Eigen::VectorXd::Zero(motifs.motifVector_.size());
            Eigen::MatrixXd interEnergies = Eigen::MatrixXd::Zero(motifs.motifVector_.size(), motifs.motifVector_.size());

            for(auto motif = motifs.motifVector_.begin(); motif != motifs.motifVector_.end(); ++motif) {
                auto motifId = std::distance(motifs.motifVector_.begin(), motif);

                // self-interaction
                intraEnergies(motifId) = calculateSelfInteractionEnergy(*motif, Te, Vee, Ven, Vnn);

                for (auto otherMotif = std::next(motif); otherMotif != motifs.motifVector_.end(); ++otherMotif) {
                    auto otherMotifId = std::distance(motifs.motifVector_.begin(), otherMotif);

                    // interaction with other motifs
                    interEnergies(motifId, otherMotifId) = caclulateInteractionEnergy(*motif, *otherMotif, Te, Vee, Ven, Vnn);
                }
            }
            return {intraEnergies, interEnergies};
        }

        double caclulateInteractionEnergy(const Motif &motif, const Motif &otherMotif,
                                          const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                          const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            double inter = 0;

            for (auto i : motif.electronIndices())
                for (auto j : otherMotif.electronIndices())
                    inter += Vee(i, j); // Vee


            for (auto i : motif.electronIndices())
                for (auto l : otherMotif.atomIndices()) {
                    inter += Ven(i, l); // Ven ->
                }

            for (auto j : otherMotif.electronIndices())
                for (auto k : motif.atomIndices())
                    inter += Ven(j, k); // Ven <-


            for (auto k : motif.atomIndices())
                for (auto l : otherMotif.atomIndices())
                    inter +=  Vnn(k, l); // Vnn

            return inter;
        }

        double calculateSelfInteractionEnergy(const Motif &motif,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn
                ) {

            double intra = 0;

            for (auto i = motif.electronIndices().begin(); i != motif.electronIndices().end(); ++i) {
                intra += Te(*i); // Te

                for (auto k = motif.atomIndices().begin(); k != motif.atomIndices().end(); ++k)
                    intra += Ven(*i, *k); // Ven

            }

            for (auto i = motif.electronIndices().begin(); i != motif.electronIndices().end(); ++i) {
                for (auto j = std::next(i); j != motif.electronIndices().end(); ++j)
                    intra += Vee(*i, *j); // Vee
            }

            // atoms in motif ( zero contribution if zero or one atom is present)
            for (auto k = motif.atomIndices().begin(); k != motif.atomIndices().end(); ++k)
                for (auto l = std::next(k); l != motif.atomIndices().end(); ++l)
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
