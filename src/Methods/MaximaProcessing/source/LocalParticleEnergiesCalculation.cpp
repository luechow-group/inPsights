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
#include <map>
#include <iterator>

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
    size_t numberOfElectrons = group.representative()->maximum().numberOfEntities();

    // add only, if the wanted selectedElectronsCount was really found within the group
    if (group.isLeaf() && size_t(group.getSelectedElectronsCount()) == selectedElectronsCount_) {

        const auto &ref = *group.representative();
        auto permutedNuclei = ref.nuclei();
        permutedNuclei.permute(ref.nuclearPermutation());
        auto VnnMat = CoulombPotential::energies(permutedNuclei);

        selectedRestInter(group, numberOfElectrons, permutedNuclei, VnnMat);
        bondEnergyCalculation(group, numberOfElectrons, permutedNuclei, VnnMat);

    } else {
        for (const auto& subgroup : group)
            add(subgroup);
    }
}

void LocalParticleEnergiesCalculator::selectedRestInter(const Group &group, size_t numberOfElectrons,
                                                        const AtomsVector &permutedNuclei,
                                                        const Eigen::MatrixXd &VnnMat) {
    // create electron indice lists
    std::vector<size_t> selectedElectronIndices, remainingElectronIndices, remainingNucleiIndices;
    createIndiceLists(numberOfElectrons, permutedNuclei,
            selectedElectronIndices,
            remainingElectronIndices,
            remainingNucleiIndices);

    for (auto &id : group.representative()->sampleIds()) {
        // add all samples
        auto &electrons = samples_[id].sample_;
        Eigen::VectorXd TeVec = samples_[id].kineticEnergies_;
        Eigen::MatrixXd VeeMat = CoulombPotential::energies(electrons);
        Eigen::MatrixXd VenMat = CoulombPotential::energies(electrons, permutedNuclei);

        // Te
        double sumTe_selected = 0.0, sumTe_intraRest = 0.0;
        for(auto idx : selectedElectronIndices)
            sumTe_selected += TeVec[idx];
        for(auto idx : remainingElectronIndices)
            sumTe_intraRest += TeVec[idx];

        localEnergies.Te.selected.add(Eigen::Matrix<double,1,1>(sumTe_selected));
        localEnergies.Te.rest.add(Eigen::Matrix<double,1,1>(sumTe_intraRest));
        localEnergies.Te.inter.add(Eigen::Matrix<double,1,1>(0.0));

        // Vee
        double sumVee_selected = 0.0, sumVee_intraRest = 0.0, sumVee_inter = 0.0;

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
                    sumVee_intraRest += VeeMat(remainingElectronIndices[i], remainingElectronIndices[j]);
        }

        localEnergies.Vee.selected.add(Eigen::Matrix<double,1,1>(sumVee_selected));
        localEnergies.Vee.rest.add(Eigen::Matrix<double,1,1>(sumVee_intraRest));
        localEnergies.Vee.inter.add(Eigen::Matrix<double,1,1>(sumVee_inter));

        // Vnn
        double sumVnn_selected = 0.0, sumVnn_intraRest = 0.0, sumVnn_inter = 0.0;

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
                    sumVnn_intraRest += VnnMat(remainingNucleiIndices[i], remainingNucleiIndices[j]);
        }

        localEnergies.Vnn.selected.add(Eigen::Matrix<double,1,1>(sumVnn_selected));
        localEnergies.Vnn.rest.add(Eigen::Matrix<double,1,1>(sumVnn_intraRest));
        localEnergies.Vnn.inter.add(Eigen::Matrix<double,1,1>(sumVnn_inter));

        // Ven
        double sumVen_selected = 0.0, sumVen_intraRest = 0.0, sumVen_inter = 0.0;
        for (auto i : selectedElectronIndices)
            for (auto k : selectedNucleiIndices_)
                sumVen_selected += VenMat(i,k);

        for (auto i : remainingElectronIndices)
            for (auto k : remainingNucleiIndices)
                sumVen_intraRest+= VenMat(i,k);

        for (auto i : selectedElectronIndices)
            for (auto k : remainingNucleiIndices)
                sumVen_inter+= VenMat(i,k);

        for (auto i : remainingElectronIndices)
            for (auto k : selectedNucleiIndices_)
                sumVen_inter+= VenMat(i,k);

        localEnergies.Ven.selected.add(Eigen::Matrix<double,1,1>(sumVen_selected));
        localEnergies.Ven.rest.add(Eigen::Matrix<double,1,1>(sumVen_intraRest));
        localEnergies.Ven.inter.add(Eigen::Matrix<double,1,1>(sumVen_inter));

        localEnergies.E.selected.add(Eigen::Matrix<double,1,1>(sumTe_selected + sumVee_selected + sumVen_selected + sumVnn_selected));
        localEnergies.E.rest.add(Eigen::Matrix<double,1,1>(sumTe_intraRest + sumVee_intraRest + sumVen_intraRest + sumVnn_intraRest));
        localEnergies.E.inter.add(Eigen::Matrix<double,1,1>(sumVee_inter + sumVen_inter + sumVnn_inter));
    }
}


void LocalParticleEnergiesCalculator::bondEnergyCalculation(
        const Group &group, size_t numberOfElectrons,
        const AtomsVector &permutedNuclei,
        const Eigen::MatrixXd &VnnMat) {

    std::vector<size_t> selectedElectronIndices, restElectronsIndices, restNucleiIndices;
    createIndiceLists(numberOfElectrons, permutedNuclei,
                      selectedElectronIndices,
                      restElectronsIndices,
                      restNucleiIndices);

    // determine core electrons
    std::vector<size_t> selectedNonCoreElectronsIndices;
    std::map<Eigen::Index, std::set<Eigen::Index>> coreElectronsMap;

    // make a map of core electrons and refine the selected electrons list
    for(auto electronIndex : selectedElectronIndices) {
        MolecularGeometry molecule(permutedNuclei, group.representative()->maximum());
        auto [atNucleusQ, nucleusIndex] = molecule.electronAtNucleusQ(electronIndex);

        if( atNucleusQ && std::find(
                std::begin(selectedNucleiIndices_),
                std::end(selectedNucleiIndices_),
                nucleusIndex) == std::end(selectedNucleiIndices_)){
            spdlog::critical("Core-electron {} found at nucleus {} which is not in the list of selected nuclei ({}).",
                    electronIndex,
                    nucleusIndex,
                    ToString::stdvectorLongUIntToString(selectedNucleiIndices_));
        }

        if(atNucleusQ &&
        (permutedNuclei[nucleusIndex].type() != Elements::ElementType::H
        || permutedNuclei[nucleusIndex].type() != Elements::ElementType::He))
             coreElectronsMap[nucleusIndex].emplace(electronIndex);
        else
            selectedNonCoreElectronsIndices.emplace_back(electronIndex);
    }

    std::vector<size_t> selectedNonCoreNucleiIndices;
    for(auto k : selectedNucleiIndices_)
        if(coreElectronsMap.find(k) == std::end(coreElectronsMap))
            selectedNonCoreNucleiIndices.emplace_back(k);
    std::sort(std::begin(selectedNonCoreNucleiIndices), std::end(selectedNonCoreNucleiIndices));

    std::sort(std::begin(selectedNonCoreElectronsIndices), std::end(selectedNonCoreElectronsIndices));



    spdlog::debug("selected non-core electrons: {}", ToString::stdvectorLongUIntToString(selectedNonCoreElectronsIndices));
    spdlog::debug("remaining electrons: {}", ToString::stdvectorLongUIntToString(restElectronsIndices));

    spdlog::debug("selected nuclei: {}", ToString::stdvectorLongUIntToString(selectedNucleiIndices_));

    for(auto [nucleusIndex, electronIndices] : coreElectronsMap) {
        std::vector<size_t>temp(electronIndices.begin(),electronIndices.end());
        spdlog::debug("\tcore: {} core electrons: {}", nucleusIndex, ToString::stdvectorLongUIntToString(temp));
    }

    spdlog::debug("selected non-core nuclei: {}", ToString::stdvectorLongUIntToString(selectedNonCoreNucleiIndices));
    spdlog::debug("remaining nuclei: {}\n", ToString::stdvectorLongUIntToString(restNucleiIndices));


    size_t counter;

    for (auto &id : group.representative()->sampleIds()) {
        // add all samples
        auto &electrons = samples_[id].sample_;
        Eigen::VectorXd TeVec = samples_[id].kineticEnergies_;
        Eigen::MatrixXd VeeMat = CoulombPotential::energies(electrons);
        Eigen::MatrixXd VenMat = CoulombPotential::energies(electrons, permutedNuclei);

        // Te
        double sumTe_intraBond = 0.0, sumTe_intraRest = 0.0;
        std::map<Eigen::Index, double> sumTe_intraCores;
        for(auto idx : selectedNonCoreElectronsIndices)
            sumTe_intraBond += TeVec[idx];

        for(auto [nucleusIndex, electronIndices] : coreElectronsMap) {
            if (sumTe_intraCores.find(nucleusIndex) == std::end(sumTe_intraCores))
                sumTe_intraCores[nucleusIndex] = 0.0;

            for(auto electronIndex : electronIndices)
                sumTe_intraCores[nucleusIndex] += TeVec[electronIndex];
        }

        for(auto idx : restElectronsIndices)
            sumTe_intraRest += TeVec[idx];

        localBondEnergies.Te.intraBond.add(Eigen::Matrix<double,1,1>(sumTe_intraBond));
        localBondEnergies.Te.intraRest.add(Eigen::Matrix<double,1,1>(sumTe_intraRest));
        localBondEnergies.Te.interBondRest.add(Eigen::Matrix<double,1,1>(0.0));

        Eigen::VectorXd sumTe_intraCoresVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumTe_interCoresBondVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumTe_interCoresRestVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::MatrixXd sumTe_interCoresCoreMat = Eigen::MatrixXd::Zero(coreElectronsMap.size(),coreElectronsMap.size());

        counter = 0;
        for (auto [nucleusIndex, electronIndices] : coreElectronsMap) {
            sumTe_intraCoresVec[counter] = sumTe_intraCores[nucleusIndex];
            counter++;
        }
        localBondEnergies.Te.intraCores.add(sumTe_intraCoresVec);
        localBondEnergies.Te.interCoresBond.add(sumTe_interCoresBondVec);
        localBondEnergies.Te.interCoresRest.add(sumTe_interCoresRestVec);
        localBondEnergies.Te.interCoresCore.add(sumTe_interCoresCoreMat);
        localBondEnergies.Te.intraBondAndInterCoresBond.add(Eigen::Matrix<double,1,1>(sumTe_intraBond));

        // Vee
        double sumVee_intraBond = 0.0, sumVee_intraRest = 0.0, sumVee_interBondRest = 0.0;
        std::map<Eigen::Index, double> sumVee_intraCores, sumVee_interCoresBond, sumVee_interCoresRest;
        std::map<Eigen::Index, std::map<Eigen::Index, double>> sumVee_interCoresCore;

        if(!selectedNonCoreElectronsIndices.empty()) {
            for (size_t i = 0; i < selectedNonCoreElectronsIndices.size() - 1; ++i)
                for (size_t j = i + 1; j < selectedNonCoreElectronsIndices.size(); ++j)
                    sumVee_intraBond += VeeMat(selectedNonCoreElectronsIndices[i], selectedNonCoreElectronsIndices[j]);
        }

        if(!restElectronsIndices.empty()) {
            for (size_t i = 0; i < restElectronsIndices.size() - 1; ++i)
                for (size_t j = i + 1; j < restElectronsIndices.size(); ++j)
                    sumVee_intraRest += VeeMat(restElectronsIndices[i], restElectronsIndices[j]);
        }

        for (auto i : selectedNonCoreElectronsIndices) {
            for (auto j : restElectronsIndices)
                sumVee_interBondRest += VeeMat(i, j);

            for(auto [nucleusIndex, electronIndices] : coreElectronsMap) {
                if (sumVee_interCoresBond.find(nucleusIndex) == std::end(sumVee_interCoresBond))
                    sumVee_interCoresBond[nucleusIndex] = 0.0;

                for(auto electronIndex : electronIndices)
                    sumVee_interCoresBond[nucleusIndex] += VeeMat(i, electronIndex);
            }
        }

        for (auto [nucleusIndex, electronIndices] : coreElectronsMap){
            if (sumVee_intraCores.find(nucleusIndex) == std::end(sumVee_intraCores))
                sumVee_intraCores[nucleusIndex] = 0.0;

            if(!electronIndices.empty()) {
                for (auto it = std::begin(electronIndices); it != std::prev(std::end(electronIndices)); it++)
                    for (auto jt = std::next(it); jt != std::end(electronIndices); jt++)
                        sumVee_intraCores[nucleusIndex] += VeeMat(*it, *jt);
            }
        }

        for (auto i : restElectronsIndices) {
            for (auto[nucleusIndex, electronIndices] : coreElectronsMap) {
                if (sumVee_interCoresRest.find(nucleusIndex) == std::end(sumVee_interCoresRest))
                    sumVee_interCoresRest[nucleusIndex] = 0.0;

                for (auto electronIndex : electronIndices)
                    sumVee_interCoresRest[nucleusIndex] += VeeMat(i, electronIndex);
            }
        }

        if(!coreElectronsMap.empty()) {
            for (auto it = std::begin(coreElectronsMap); it != std::prev(std::end(coreElectronsMap)); it++) {
                for (auto jt = std::next(it); jt != std::end(coreElectronsMap); jt++) {

                    if(sumVee_interCoresCore.find(it->first) == std::end(sumVee_interCoresCore))
                        sumVee_interCoresCore[it->first] = {};

                    if(sumVee_interCoresCore[it->first].find(jt->first) == std::end(sumVee_interCoresCore[it->first]))
                        sumVee_interCoresCore[it->first][jt->first] = 0.0;

                    for (auto i : it->second)
                        for (auto j : jt->second) {
                            sumVee_interCoresCore[it->first][jt->first] += VeeMat(i, j);
                        }
                }
            }
        }
        localBondEnergies.Vee.intraBond.add(Eigen::Matrix<double,1,1>(sumVee_intraBond));
        localBondEnergies.Vee.intraRest.add(Eigen::Matrix<double,1,1>(sumVee_intraRest));
        localBondEnergies.Vee.interBondRest.add(Eigen::Matrix<double,1,1>(sumVee_interBondRest));

        Eigen::VectorXd sumVee_intraCoresVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumVee_interCoresBondVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumVee_interCoresRestVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::MatrixXd sumVee_interCoresCoreMat = Eigen::MatrixXd::Zero(coreElectronsMap.size(),coreElectronsMap.size());


        double sumVee_intraBondAndInterCoresBond = sumVee_intraBond;
        counter = 0;
        for (auto [nucleusIndex, electronIndices] : coreElectronsMap) {

            sumVee_intraCoresVec[counter] = sumVee_intraCores[nucleusIndex];
            sumVee_interCoresBondVec[counter] = sumVee_interCoresBond[nucleusIndex];
            sumVee_intraBondAndInterCoresBond += sumVee_interCoresBond[nucleusIndex];
            sumVee_interCoresRestVec[counter] = sumVee_interCoresRest[nucleusIndex];

            size_t counter2 = 0;
            for (auto [nucleusIndex2, electronIndices2] : coreElectronsMap) {
                sumVee_interCoresCoreMat(counter,counter2) = sumVee_interCoresCore[nucleusIndex][nucleusIndex2];
                counter2++;
            }
            counter++;
        }
        localBondEnergies.Vee.intraCores.add(sumVee_intraCoresVec);
        localBondEnergies.Vee.interCoresBond.add(sumVee_interCoresBondVec);
        localBondEnergies.Vee.interCoresRest.add(sumVee_interCoresRestVec);
        localBondEnergies.Vee.interCoresCore.add(sumVee_interCoresCoreMat);
        localBondEnergies.Vee.intraBondAndInterCoresBond.add(Eigen::Matrix<double,1,1>(sumVee_intraBondAndInterCoresBond));


        // Vnn
        double sumVnn_intraBond = 0.0, sumVnn_intraRest = 0.0, sumVnn_interBondRest = 0.0;
        std::map<Eigen::Index, double> sumVnn_interCoresBond, sumVnn_interCoresRest;
        std::map<Eigen::Index, std::map<Eigen::Index, double>> sumVnn_interCoresCore;

        if(!selectedNonCoreNucleiIndices.empty()) {
            for (size_t i = 0; i < selectedNonCoreNucleiIndices.size() - 1; ++i)
                for (size_t j = i + 1; j < selectedNonCoreNucleiIndices.size(); ++j)
                    sumVnn_intraBond += VnnMat(selectedNonCoreNucleiIndices[i], selectedNonCoreNucleiIndices[j]);
        }
        if(!restNucleiIndices.empty()) {
            for (size_t i = 0; i < restNucleiIndices.size() - 1; ++i)
                for (size_t j = i + 1; j < restNucleiIndices.size(); ++j)
                    sumVnn_intraRest += VnnMat(restNucleiIndices[i], restNucleiIndices[j]);
        }

        for (auto k : selectedNonCoreNucleiIndices)
            for (auto l : restNucleiIndices)
                sumVnn_interBondRest += VnnMat(k, l);


        for (auto k : selectedNonCoreNucleiIndices) {
            for(auto [nucleusIndex, electronIndices] : coreElectronsMap) {
                assert(long(k) != nucleusIndex);
                if (sumVnn_interCoresBond.find(nucleusIndex) == std::end(sumVnn_interCoresBond))
                    sumVnn_interCoresBond[nucleusIndex] = 0.0;

                sumVnn_interCoresBond[nucleusIndex] += VnnMat(k, nucleusIndex);
            }
        }

        for (auto [coreNucleusIndex, electronIndices] : coreElectronsMap) {
            for (auto l : restNucleiIndices) {

                if(sumVnn_interCoresRest.find(coreNucleusIndex) == std::end(sumVnn_interCoresRest))
                    sumVnn_interCoresRest[coreNucleusIndex] = 0.0;

                sumVnn_interCoresRest[coreNucleusIndex] += VnnMat(coreNucleusIndex, l);
            }
        }

        if(!coreElectronsMap.empty()) {
            for (auto it = std::begin(coreElectronsMap); it != std::prev(std::end(coreElectronsMap)); it++) {
                for (auto jt = std::next(it); jt != std::end(coreElectronsMap); jt++) {

                    if(sumVnn_interCoresCore.find(it->first) == std::end(sumVnn_interCoresCore))
                        sumVnn_interCoresCore[it->first] = {};

                    if(sumVnn_interCoresCore[it->first].find(jt->first) == std::end(sumVnn_interCoresCore[it->first]))
                        sumVnn_interCoresCore[it->first][jt->first] = 0.0;

                    sumVnn_interCoresCore[it->first][jt->first] += VnnMat(it->first, jt->first);
                }
            }
        }

        localBondEnergies.Vnn.intraBond.add(Eigen::Matrix<double,1,1>(sumVnn_intraBond));
        localBondEnergies.Vnn.intraRest.add(Eigen::Matrix<double,1,1>(sumVnn_intraRest));
        localBondEnergies.Vnn.interBondRest.add(Eigen::Matrix<double,1,1>(sumVnn_interBondRest));

        Eigen::VectorXd sumVnn_interCoresBondVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumVnn_interCoresRestVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::MatrixXd sumVnn_interCoresCoreMat = Eigen::MatrixXd::Zero(coreElectronsMap.size(), coreElectronsMap.size());

        counter = 0;
        double sumVnn_intraBondAndInterCoresBond = sumVnn_intraBond;
        for (auto [nucleusIndex, electronIndices] : coreElectronsMap) {

            sumVnn_interCoresBondVec[counter] = sumVnn_interCoresBond[nucleusIndex];
            sumVnn_intraBondAndInterCoresBond += sumVnn_interCoresBond[nucleusIndex];
            sumVnn_interCoresRestVec[counter] = sumVnn_interCoresRest[nucleusIndex];

            size_t counter2 = 0;
            for (auto [nucleusIndex2, electronIndice2] : coreElectronsMap) {
                sumVnn_interCoresCoreMat(counter,counter2) = sumVnn_interCoresCore[nucleusIndex][nucleusIndex2];
                counter2++;
            }
            counter++;
        }

        localBondEnergies.Vnn.interCoresBond.add(sumVnn_interCoresBondVec);
        localBondEnergies.Vnn.interCoresRest.add(sumVnn_interCoresRestVec);
        localBondEnergies.Vnn.interCoresCore.add(sumVnn_interCoresCoreMat);
        localBondEnergies.Vnn.intraBondAndInterCoresBond.add(Eigen::Matrix<double,1,1>(sumVnn_intraBondAndInterCoresBond));

        // Ven
        double sumVen_intraBond = 0.0, sumVen_intraRest = 0.0, sumVen_interBondRest = 0.0;
        std::map<Eigen::Index, double> sumVen_intraCores, sumVen_interCoresBond, sumVen_interCoresRest;
        std::map<Eigen::Index, std::map<Eigen::Index, double>> sumVen_interCoresCore;

        for (auto i : selectedNonCoreElectronsIndices)
            for (auto k : selectedNonCoreNucleiIndices)
                sumVen_intraBond += VenMat(i, k);

        for (auto restElectronIndex : restElectronsIndices)
            for (auto restNucleusIndex : restNucleiIndices)
                sumVen_intraRest += VenMat(restElectronIndex, restNucleusIndex);

        for (auto i : selectedNonCoreElectronsIndices)
            for (auto k : restNucleiIndices)
                sumVen_interBondRest += VenMat(i, k);
        for (auto i : restElectronsIndices)
            for (auto k : selectedNonCoreNucleiIndices)
                sumVen_interBondRest += VenMat(i, k);

        for (auto [nucleusIndex, electronIndices] : coreElectronsMap){
            sumVen_intraCores[nucleusIndex] = 0.0;

            for(auto electronIndex : electronIndices)
                sumVen_intraCores[nucleusIndex] += VenMat(electronIndex, nucleusIndex);
        }

        for(auto [nucleusIndex, electronIndices] : coreElectronsMap) {
            sumVen_interCoresBond[nucleusIndex] = 0.0;

            // core-nuclei of selection with selected non-core electrons
            for (auto i : selectedNonCoreElectronsIndices)
                sumVen_interCoresBond[nucleusIndex] += VenMat(i, nucleusIndex);

            // selected non-core nuclei with core-electrons
            for (auto i : electronIndices)
                for(auto k : selectedNonCoreNucleiIndices)
                sumVen_interCoresBond[nucleusIndex] += VenMat(i, k);
        }

        for (auto[nucleusIndex, coreElectronIndices] : coreElectronsMap) {
            if (sumVen_interCoresRest.find(nucleusIndex) == std::end(sumVen_interCoresRest))
                sumVen_interCoresRest[nucleusIndex] = 0.0;

            for (auto i : coreElectronIndices)
                for(auto k : restNucleiIndices)
                sumVen_interCoresRest[nucleusIndex] += VenMat(i, k);

            for (auto i : restElectronsIndices)
                sumVen_interCoresRest[nucleusIndex] += VenMat(i, nucleusIndex);
        }

        if(!coreElectronsMap.empty()) {
            for (auto kt = std::begin(coreElectronsMap); kt != std::prev(std::end(coreElectronsMap)); kt++) {
                for (auto lt = std::next(kt); lt != std::end(coreElectronsMap); lt++) {

                    if(sumVen_interCoresCore.find(kt->first) == std::end(sumVen_interCoresCore))
                        sumVen_interCoresCore[kt->first] = {};

                    if(sumVen_interCoresCore[kt->first].find(lt->first) == std::end(sumVen_interCoresCore[kt->first]))
                        sumVen_interCoresCore[kt->first][lt->first] = 0.0;


                    // interaction of electrons i of core k with core l
                    for (auto i : kt->second)
                        sumVen_interCoresCore[kt->first][lt->first] += VenMat(i, lt->first);

                    // interaction of electrons j of core l with core k
                    for (auto j : lt->second)
                        sumVen_interCoresCore[kt->first][lt->first] += VenMat(j, kt->first);
                }
            }
        }


        localBondEnergies.Ven.intraBond.add(Eigen::Matrix<double,1,1>(sumVen_intraBond));
        localBondEnergies.Ven.intraRest.add(Eigen::Matrix<double,1,1>(sumVen_intraRest));
        localBondEnergies.Ven.interBondRest.add(Eigen::Matrix<double,1,1>(sumVen_interBondRest));

        Eigen::VectorXd sumVen_intraCoresVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumVen_interCoresBondVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumVen_interCoresRestVec = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::MatrixXd sumVen_interCoresCoreMat = Eigen::MatrixXd::Zero(coreElectronsMap.size(),coreElectronsMap.size());

        counter = 0;
        double sumVen_intraBondAndInterCoresBond = sumVen_intraBond;
        for (auto [nucleusIndex, electronIndices] : coreElectronsMap) {

            sumVen_intraCoresVec[counter] = sumVen_intraCores[nucleusIndex];
            sumVen_interCoresBondVec[counter] = sumVen_interCoresBond[nucleusIndex];
            sumVen_intraBondAndInterCoresBond += sumVen_interCoresBond[nucleusIndex];
            sumVen_interCoresRestVec[counter] = sumVen_interCoresRest[nucleusIndex];

            size_t counter2 = 0;
            for (auto [nucleusIndex2, electronIndices2] : coreElectronsMap) {
                sumVen_interCoresCoreMat(counter,counter2) = sumVen_interCoresCore[nucleusIndex][nucleusIndex2];
                counter2++;
            }
            counter++;
        }

        localBondEnergies.Ven.intraCores.add(sumVen_intraCoresVec);
        localBondEnergies.Ven.interCoresBond.add(sumVen_interCoresBondVec);
        localBondEnergies.Ven.interCoresRest.add(sumVen_interCoresRestVec);
        localBondEnergies.Ven.interCoresCore.add(sumVen_interCoresCoreMat);
        localBondEnergies.Ven.intraBondAndInterCoresBond.add(Eigen::Matrix<double,1,1>(sumVen_intraBondAndInterCoresBond));


        // E
        localBondEnergies.E.intraBond.add(Eigen::Matrix<double,1,1>(sumTe_intraBond + sumVee_intraBond + sumVen_intraBond + sumVnn_intraBond ));
        localBondEnergies.E.intraRest.add(Eigen::Matrix<double,1,1>(sumTe_intraRest + sumVee_intraRest + sumVen_intraRest + sumVnn_intraRest));
        localBondEnergies.E.interBondRest.add(Eigen::Matrix<double,1,1>(sumVee_interBondRest + sumVen_interBondRest + sumVnn_interBondRest));

        Eigen::VectorXd sumE_intraCores = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumE_interCoresBond = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::VectorXd sumE_interCoresRest = Eigen::VectorXd::Zero(coreElectronsMap.size());
        Eigen::MatrixXd sumE_interCoresCore = Eigen::MatrixXd::Zero(coreElectronsMap.size(),coreElectronsMap.size());

        counter = 0;
        for (auto [k, electronIndices] : coreElectronsMap) {

            sumE_intraCores[counter] = sumTe_intraCores[k] + sumVee_intraCores[k] + sumVen_intraCores[k];
            sumE_interCoresBond[counter] = sumVee_interCoresBond[k] + sumVen_interCoresBond[k] + sumVnn_interCoresBond[k];
            sumE_interCoresRest[counter] = sumVee_interCoresRest[k] + sumVen_interCoresRest[k] + sumVnn_interCoresRest[k];


            size_t counter2 = 0;
            for (auto [l, electronIndice2] : coreElectronsMap) {
                sumE_interCoresCore(counter,counter2) =
                        sumVee_interCoresCore[k][l]+ sumVen_interCoresCore[k][l]+ sumVnn_interCoresCore[k][l];

                counter2++;
            }
            counter++;
        }
        localBondEnergies.E.intraBondAndInterCoresBond.add(Eigen::Matrix<double,1,1>(
                sumTe_intraBond +
                sumVee_intraBondAndInterCoresBond +
                sumVen_intraBondAndInterCoresBond +
                sumVnn_intraBondAndInterCoresBond));
        localBondEnergies.E.intraCores.add(sumE_intraCores);
        localBondEnergies.E.interCoresBond.add(sumE_interCoresBond);
        localBondEnergies.E.interCoresRest.add(sumE_interCoresRest);
        localBondEnergies.E.interCoresCore.add(sumE_interCoresCore);
    }
}

void LocalParticleEnergiesCalculator::createIndiceLists(size_t numberOfElectrons, const AtomsVector &permutedNuclei,
                                                        std::vector<size_t> &selectedElectronIndices,
                                                        std::vector<size_t> &remainingElectronIndices,
                                                        std::vector<size_t> &remainingNucleiIndices) const {// create electron indice lists
    for (size_t i = 0; i < selectedElectronsCount_; ++i)
        selectedElectronIndices.emplace_back(i);
    for (size_t i = selectedElectronsCount_; i < numberOfElectrons; ++i)
        remainingElectronIndices.emplace_back(i);

    // create remaining nuclei indice list
    for (long k = 0; k < permutedNuclei.numberOfEntities(); ++k)
        if(std::find(std::begin(selectedNucleiIndices_), std::end(selectedNucleiIndices_), k) == std::end(
                selectedNucleiIndices_))
            remainingNucleiIndices.emplace_back(k);
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

    Emitter &operator<<(Emitter &out, const LocalParticleEnergiesCalculator::LocalBondEnergyResults &rhs) {
        out << BeginMap
            << Key << "bondEnergy" << Value << rhs.intraBondAndInterCoresBond << Newline
            << Key << "intraBond" << Value << rhs.intraBond << Newline
            << Key << "intraRest" << Value << rhs.intraRest << Newline
            << Key << "interBondRest" << Value << rhs.interBondRest << Newline
            << Key << "intraCores" << Value << rhs.intraCores << Newline
            << Key << "interCoresBond" << Value << rhs.interCoresBond << Newline
            << Key << "interCoresRest" << Value << rhs.interCoresRest << Newline
            << Key << "interCoresCore" << Value << rhs.interCoresCore << Newline
            << EndMap;
        return out;
    }


    Emitter &operator<<(Emitter &out, const LocalParticleEnergiesCalculator &rhs) {
        out << BeginMap
            << Key << "LocalEnergies" << Value << rhs.localEnergies
             << Key << "LocalBondEnergies" << Value << rhs.localBondEnergies
             << EndMap;
        return out;
    }
}