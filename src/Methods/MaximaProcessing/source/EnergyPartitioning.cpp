/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <EnergyPartitioning.h>
#include <algorithm>
#include <Motif.h>
#include <Enumerate.h>

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

    namespace MolecularSelectionBased{

        std::pair<Eigen::VectorXd, Eigen::MatrixXd> calculateInteractionEnergies(const Motifs &motifs,
                                                                                 const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                                                 const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            std::vector<MolecularSelection> selections(motifs.motifs_.size());
            for(auto [i, m] : enumerate(motifs.motifs_))
                selections[i] = static_cast<MolecularSelection>(m);

            return calculateInteractionEnergies(selections, Te, Vee, Ven, Vnn);
        }

        std::pair<Eigen::VectorXd, Eigen::MatrixXd> calculateInteractionEnergies(const std::vector<MolecularSelection> &selections,
                                                                                 const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                                                 const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {

            Eigen::VectorXd intraEnergies = Eigen::VectorXd::Zero(selections.size());
            Eigen::MatrixXd interEnergies = Eigen::MatrixXd::Zero(selections.size(), selections.size());

            for(auto motif = selections.begin(); motif != selections.end(); ++motif) {
                auto motifId = std::distance(selections.begin(), motif);

                // self-interaction
                intraEnergies(motifId) = calculateSelfInteractionEnergy(*motif, Te, Vee, Ven, Vnn);

                for (auto otherMotif = std::next(motif); otherMotif != selections.end(); ++otherMotif) {
                    auto otherMotifId = std::distance(selections.begin(), otherMotif);

                    // interaction with other selections
                    interEnergies(motifId, otherMotifId) = caclulateInteractionEnergy(*motif, *otherMotif, Te, Vee, Ven, Vnn);
                }
            }

            // symmetrize the matrix of inter energies
            interEnergies = interEnergies.selfadjointView<Eigen::Upper>();

            return {intraEnergies, interEnergies};
        }

        double caclulateInteractionEnergy(const MolecularSelection &selA, const MolecularSelection &selB,
                                          const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                          const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            double inter = 0;

            for (auto i : selA.electrons_.indices())
                for (auto j : selB.electrons_.indices())
                    inter += Vee(i, j); // Vee


            for (auto i : selA.electrons_.indices())
                for (auto l : selB.nuclei_.indices()) {
                    inter += Ven(i, l); // Ven ->
                }

            for (auto j : selB.electrons_.indices())
                for (auto k : selA.nuclei_.indices())
                    inter += Ven(j, k); // Ven <-


            for (auto k : selA.nuclei_.indices())
                for (auto l : selB.nuclei_.indices())
                    inter +=  Vnn(k, l); // Vnn

            return inter;
        }

        double calculateSelfInteractionEnergy(const MolecularSelection &sel,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn
                ) {

            double intra = 0;

            for (auto i = sel.electrons_.indices().begin(); i != sel.electrons_.indices().end(); ++i) {
                intra += Te(*i); // Te

                for (auto k = sel.nuclei_.indices().begin(); k != sel.nuclei_.indices().end(); ++k)
                    intra += Ven(*i, *k); // Ven

                for (auto j = std::next(i); j != sel.electrons_.indices().end(); ++j)
                    intra += Vee(*i, *j); // Vee
            }

            // atoms in motif ( zero contribution if zero or one atom is present)
            for (auto k = sel.nuclei_.indices().begin(); k != sel.nuclei_.indices().end(); ++k)
                for (auto l = std::next(k); l != sel.nuclei_.indices().end(); ++l)
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
