/* Copyright (C) 2020 Michael Heuer.
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

#include "LocalParticleEnergiesCalculation.h"
#include <CoulombPotential.h>
#include <Reference.h>

LocalParticleEnergiesCalculator::LocalParticleEnergiesCalculator(
        const std::vector<Sample> &samples,
        const AtomsVector &atoms,
        const std::vector<size_t> &nucleiIndices,
        size_t selectedElectronsCount)
        :
        samples_(samples),
        selectedNucleiIndices_(nucleiIndices),
        selectedElectronsCount_(selectedElectronsCount)
        {}

void LocalParticleEnergiesCalculator::add(const Group &group) {
    auto numberOfElectrons = group.representative()->maximum().numberOfEntities();

    if (group.isLeaf() && group.getSelectedElectronsCount() == size_t(selectedElectronsCount_)) {

        const auto &ref = *group.representative();
        auto permutedNuclei = ref.nuclei();
        permutedNuclei.permute(ref.nuclearPermutation());
        auto VnnMat = CoulombPotential::energies(permutedNuclei);


        // create electron indice lists
        std::vector<size_t> selectedElectronIndices, remainingElectronIndices, remainingNucleiIndices;
        for (size_t i = 0; i < selectedElectronsCount_; ++i)
            selectedElectronIndices.emplace_back(i);
        for (size_t i = selectedElectronsCount_; i < numberOfElectrons; ++i)
            remainingElectronIndices.emplace_back(i);

        // create remaining nuclei indice list
        for (long k = 0; k < permutedNuclei.numberOfEntities(); ++k)
            if(std::find(std::begin(selectedNucleiIndices_), std::end(selectedNucleiIndices_), k) == std::end(selectedNucleiIndices_))
                remainingNucleiIndices.emplace_back(k);


        for (auto &id : group.representative()->sampleIds()) {
            // add all samples
            auto &electrons = samples_[id].sample_;
            Eigen::VectorXd TeVec = samples_[id].kineticEnergies_;
            Eigen::MatrixXd VeeMat = CoulombPotential::energies(electrons);
            Eigen::MatrixXd VenMat = CoulombPotential::energies(electrons, permutedNuclei);

            // Te
            double sumTe_selected = 0.0, sumTe_rest = 0.0;
            for(auto idx : selectedElectronIndices)
                sumTe_selected += TeVec[idx];
            for(auto idx : remainingElectronIndices)
                sumTe_rest += TeVec[idx];

            Te.selected.add(Eigen::Matrix<double,1,1>(sumTe_selected));
            Te.rest.add(Eigen::Matrix<double,1,1>(sumTe_rest));
            Te.inter.add(Eigen::Matrix<double,1,1>(0.0));

            // Vee
            double sumVee_selected = 0.0, sumVee_rest = 0.0, sumVee_inter = 0.0;

            if(!selectedElectronIndices.empty()) {
                for (size_t i = 0; i < selectedElectronIndices.size() - 1; ++i)
                    for (size_t j = i + 1; j < selectedElectronIndices.size(); ++j)
                        sumVee_selected += VeeMat(selectedElectronIndices[i], selectedElectronIndices[j]);
            }

            for (auto i : selectedElectronIndices)
                for (auto j : remainingElectronIndices)
                    sumVee_inter += VeeMat(i, j);

            if(!remainingElectronIndices.empty()) {
                for (size_t i = 0; i < remainingElectronIndices.size() - 1; ++i)
                    for (size_t j = i + 1; j < remainingElectronIndices.size(); ++j)
                        sumVee_rest += VeeMat(remainingElectronIndices[i], remainingElectronIndices[j]);
            }

            Vee.selected.add(Eigen::Matrix<double,1,1>(sumVee_selected));
            Vee.rest.add(Eigen::Matrix<double,1,1>(sumVee_rest));
            Vee.inter.add(Eigen::Matrix<double,1,1>(sumVee_inter));

            // Vnn
            double sumVnn_selected = 0.0, sumVnn_rest = 0.0, sumVnn_inter = 0.0;

            if(!selectedNucleiIndices_.empty()) {
                for (size_t i = 0; i < selectedNucleiIndices_.size() - 1; ++i)
                    for (size_t j = i + 1; j < selectedNucleiIndices_.size(); ++j)
                        sumVnn_selected += VnnMat(selectedNucleiIndices_[i], selectedNucleiIndices_[j]);
            }

            for (auto i : selectedNucleiIndices_)
                for (auto j : remainingNucleiIndices)
                    sumVnn_inter += VnnMat(i, j);

            if(!remainingNucleiIndices.empty()) {
                for (size_t i = 0; i < remainingNucleiIndices.size() - 1; ++i)
                    for (size_t j = i + 1; j < remainingNucleiIndices.size(); ++j)
                        sumVnn_rest += VnnMat(remainingNucleiIndices[i], remainingNucleiIndices[j]);
            }

            Vnn.selected.add(Eigen::Matrix<double,1,1>(sumVnn_selected));
            Vnn.rest.add(Eigen::Matrix<double,1,1>(sumVnn_rest));
            Vnn.inter.add(Eigen::Matrix<double,1,1>(sumVnn_inter));

            // Ven
            double sumVen_selected = 0.0, sumVen_rest = 0.0, sumVen_inter = 0.0;
            for (auto i : selectedElectronIndices)
                for (auto k : selectedNucleiIndices_)
                    sumVen_selected += VenMat(i,k);

            for (auto i : remainingElectronIndices)
                for (auto k : remainingNucleiIndices)
                    sumVen_rest+= VenMat(i,k);

            for (auto i : selectedElectronIndices)
                for (auto k : remainingNucleiIndices)
                    sumVen_inter+= VenMat(i,k);

            for (auto i : remainingElectronIndices)
                for (auto k : selectedNucleiIndices_)
                    sumVen_inter+= VenMat(i,k);

            Ven.selected.add(Eigen::Matrix<double,1,1>(sumVen_selected));
            Ven.rest.add(Eigen::Matrix<double,1,1>(sumVen_rest));
            Ven.inter.add(Eigen::Matrix<double,1,1>(sumVen_inter));

            E.selected.add(Eigen::Matrix<double,1,1>(
                    sumTe_selected + sumVee_selected + sumVen_selected + sumVnn_selected));
            E.rest.add(Eigen::Matrix<double,1,1>(
                    sumTe_rest + sumVee_rest + sumVen_rest + sumVnn_rest));
            E.inter.add(Eigen::Matrix<double,1,1>(sumVee_inter + sumVen_inter + sumVnn_inter));
        }
    } else {
        for (const auto& subgroup : group)
            add(subgroup);
    }
}

namespace YAML {

    Emitter &operator<<(Emitter &out, const LocalParticleEnergiesCalculator::LocalEnergyResults &rhs) {
        out << BeginMap
            << Key << "selected" << Value << rhs.selected << Newline
            << Key << "rest" << Value << rhs.rest << Newline
            << Key << "inter" << Value << rhs.inter
            << EndMap;
        return out;
    }

    Emitter &operator<<(Emitter &out, const LocalParticleEnergiesCalculator &rhs) {
        out << BeginMap
            << Key << "E" << Value << rhs.E << Newline
            << Key << "Te" << Value << rhs.Te << Newline
            << Key << "Vee" << Value << rhs.Vee << Newline
            << Key << "Ven" << Value << rhs.Ven << Newline
            << Key << "Vnn" << Value << rhs.Vnn
            << EndMap;
        return out;
    }
}