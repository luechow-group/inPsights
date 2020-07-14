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

#ifndef INPSIGHTS_LOCALPARTICLEENERGIESCALCULATION_H
#define INPSIGHTS_LOCALPARTICLEENERGIESCALCULATION_H

#include <Statistics.h>
#include <Sample.h>
#include <Cluster.h>
#include <ParticleSelection.h>

class LocalParticleEnergiesCalculator {
public:
    struct LocalEnergyResults {
        SingleValueStatistics selected, rest, inter;
    };

    struct LocalBondEnergyResults {
        SingleValueStatistics intraBondAndInterCoresBond,intraBond, intraRest, interBondRest;
        VectorStatistics intraCores, interCoresBond, interCoresRest;
        TriangularMatrixStatistics interCoresCore;
    };

    template<typename Type>
    struct ResultsBundle {
        Type E, Te, Vee, Ven, Vnn;
    };

    LocalParticleEnergiesCalculator(
            const std::vector<Sample> &samples,
            const AtomsVector &atoms,
            const std::vector<size_t> &nucleiIndices,
            size_t selectedElectronsCount);

    void add(const Cluster &group);

    const std::vector<Sample> &samples_;
    std::vector<size_t> selectedNucleiIndices_;
    size_t selectedElectronsCount_;

    ResultsBundle<LocalEnergyResults> localEnergies;
    ResultsBundle<LocalBondEnergyResults> localBondEnergies;


    void selectedRestInter(const Cluster &group, size_t numberOfElectrons, const AtomsVector &permutedNuclei,
            const Eigen::MatrixXd &VnnMat);

    void bondEnergyCalculation(const Cluster &group, size_t numberOfElectrons, const AtomsVector &permutedNuclei,
                               const Eigen::MatrixXd &VnnMat) ;

    void createIndiceLists(size_t numberOfElectrons, const AtomsVector &permutedNuclei,
                           std::vector<size_t> &selectedElectronIndices, std::vector<size_t> &remainingElectronIndices,
                           std::vector<size_t> &remainingNucleiIndices) const;
};
namespace YAML {
    class Emitter;

    Emitter& operator<< (Emitter& out, const LocalParticleEnergiesCalculator::LocalEnergyResults& rhs);
    Emitter& operator<< (Emitter& out, const LocalParticleEnergiesCalculator::LocalBondEnergyResults& rhs);


    template<typename Type>
    Emitter &operator<<(Emitter &out, const LocalParticleEnergiesCalculator::ResultsBundle<Type> &rhs) {
        out << BeginMap
            << Key << "E" << Comment("[Eh]") << Value << rhs.E
            << Key << "Te" << Comment("[Eh]") << Value << rhs.Te
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.Vee
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.Ven
            << Key << "Vnn" << Comment("[Eh]") << Value << rhs.Vnn
            << EndMap;
        return out;
    }

    Emitter& operator<< (Emitter& out, const LocalParticleEnergiesCalculator& rhs);


}


#endif //INPSIGHTS_LOCALPARTICLEENERGIESCALCULATION_H
