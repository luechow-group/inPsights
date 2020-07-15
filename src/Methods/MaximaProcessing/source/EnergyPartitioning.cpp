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

        std::pair<EnergyResultsBundle<Eigen::VectorXd>, EnergyResultsBundle<Eigen::MatrixXd>> calculateInteractionEnergies(const Motifs &motifs,
                                                                                 const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                                                 const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            std::vector<MolecularSelection> selections(motifs.motifs_.size());
            for(auto [i, m] : enumerate(motifs.motifs_))
                selections[i] = static_cast<MolecularSelection>(m);

            return calculateInteractionEnergies(selections, Te, Vee, Ven, Vnn);
        }

        std::pair<EnergyResultsBundle<Eigen::VectorXd>, EnergyResultsBundle<Eigen::MatrixXd>> calculateInteractionEnergies(const std::vector<MolecularSelection> &selections,
                                                                                 const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                                                                                 const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            EnergyResultsBundle<Eigen::VectorXd> intra;
            EnergyResultsBundle<Eigen::MatrixXd> inter;

            intra.init(Eigen::VectorXd::Zero(selections.size()));
            inter.init(Eigen::MatrixXd::Zero(selections.size(), selections.size()));


            for(auto motif = selections.begin(); motif != selections.end(); ++motif) {
                auto motifId = std::distance(selections.begin(), motif);

                // self-interaction
                auto intraValues = calculateSelfInteractionEnergyBundle(*motif, Te, Vee, Ven, Vnn);
                intra.E(motifId) = intraValues.E;
                intra.Te(motifId) = intraValues.Te;
                intra.Vee(motifId) = intraValues.Vee;
                intra.Ven(motifId) = intraValues.Ven;
                intra.Vnn(motifId) = intraValues.Vnn;


                for (auto otherMotif = std::next(motif); otherMotif != selections.end(); ++otherMotif) {
                    auto otherMotifId = std::distance(selections.begin(), otherMotif);

                    // interaction with other selections
                    auto interValues = calculateInteractionEnergyBundle(*motif, *otherMotif, Vee, Ven, Vnn);
                    inter.E(motifId, otherMotifId) = interValues.E;
                    inter.Vee(motifId, otherMotifId) = interValues.Vee;
                    inter.Ven(motifId, otherMotifId) = interValues.Ven;
                    inter.Vnn(motifId, otherMotifId) = interValues.Vnn;

                }
            }

            // symmetrize the matrix of inter energies
            inter.E = inter.E.selfadjointView<Eigen::Upper>();
            inter.Vee = inter.Vee.selfadjointView<Eigen::Upper>();
            inter.Ven = inter.Ven.selfadjointView<Eigen::Upper>();
            inter.Vnn = inter.Vnn.selfadjointView<Eigen::Upper>();

            return {intra, inter};
        }

        double calculateIntraTe(const MolecularSelection &sel, const Eigen::VectorXd &Te) {
            double intra = 0;

            for (auto i = sel.electrons_.indices().begin(); i != sel.electrons_.indices().end(); ++i)
                intra += Te(*i);

            return intra;
        }

        double calculateIntraVee(const MolecularSelection &sel, const Eigen::MatrixXd &Vee) {
            double intra = 0;

            for (auto i = sel.electrons_.indices().begin(); i != sel.electrons_.indices().end(); ++i)
                for (auto j = std::next(i); j != sel.electrons_.indices().end(); ++j)
                    intra += Vee(*i, *j);

            return intra;
        }

        double
        calculateInterVee(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Vee) {
            double inter = 0;

            for (auto i : selA.electrons_.indices())
                for (auto j : selB.electrons_.indices())
                    inter += Vee(i, j);

            return inter;
        }

        double calculateIntraVen(const MolecularSelection &sel, const Eigen::MatrixXd &Ven) {
            double intra = 0;

            for (auto i = sel.electrons_.indices().begin(); i != sel.electrons_.indices().end(); ++i)
                for (auto k = sel.nuclei_.indices().begin(); k != sel.nuclei_.indices().end(); ++k)
                    intra += Ven(*i, *k);

            return intra;
        }

        double
        calculateInterVen(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Ven) {
            double inter = 0;

            for (auto i : selA.electrons_.indices())
                for (auto l : selB.nuclei_.indices())
                    inter += Ven(i, l); // A -> B

            for (auto j : selB.electrons_.indices())
                for (auto k : selA.nuclei_.indices())
                    inter += Ven(j, k); // A <- B

            return inter;
        }

        double calculateIntraVnn(const MolecularSelection &sel, const Eigen::MatrixXd &Vnn) {
            double intra = 0;
            for (auto k = sel.nuclei_.indices().begin(); k != sel.nuclei_.indices().end(); ++k)
                for (auto l = std::next(k); l != sel.nuclei_.indices().end(); ++l)
                    intra += Vnn(*k, *l);

            return intra;
        }

        double
        calculateInterVnn(const MolecularSelection &selA, const MolecularSelection &selB, const Eigen::MatrixXd &Vnn) {
            double inter = 0;
            for (auto k : selA.nuclei_.indices())
                for (auto l : selB.nuclei_.indices())
                    inter +=  Vnn(k, l);

            return inter;
        }

        EnergyResultsBundle<double> calculateInteractionEnergyBundle(const MolecularSelection &selA, const MolecularSelection &selB,
                                          const Eigen::MatrixXd &Vee, const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn) {
            auto
            vee = calculateInterVee(selA, selB, Vee),
            ven = calculateInterVen(selA, selB, Ven),
            vnn = calculateInterVnn(selA, selB, Vnn);

            return {vee + ven + vnn, 0, vee, ven, vnn};
        }

        EnergyResultsBundle<double> calculateSelfInteractionEnergyBundle(const MolecularSelection &sel,
                const Eigen::VectorXd &Te, const Eigen::MatrixXd &Vee,
                const Eigen::MatrixXd &Ven, const Eigen::MatrixXd &Vnn
                ) {
            auto
            te = calculateIntraTe(sel, Te),
            vee = calculateIntraVee(sel, Vee),
            ven = calculateIntraVen(sel, Ven),
            vnn = calculateIntraVnn(sel, Vnn);

            return {te + vee + ven + vnn, te, vee, ven, vnn};
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
