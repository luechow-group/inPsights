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
#include <Group.h>
#include <ParticleSelection.h>

class LocalParticleEnergiesCalculator {
public:
    struct LocalEnergyResults {
        SingleValueStatistics selected, rest, inter;
    };

    LocalParticleEnergiesCalculator(
            const std::vector<Sample> &samples,
            const AtomsVector &atoms,
            const std::vector<size_t> &nucleiIndices,
            size_t selectedElectronsCount);

    void add(const Group &group);

    const std::vector<Sample> &samples_;
    std::vector<size_t> selectedNucleiIndices_;
    size_t selectedElectronsCount_;

    LocalEnergyResults E, Te, Vee, Ven, Vnn;

};
namespace YAML {
    class Node;
    class Emitter;

    template<typename Type>
    struct convert;

    //template<>
    //struct convert<LocalParticleEnergiesCalculator::LocalEnergyResults> {
    //    static Node encode(const LocalParticleEnergiesCalculator::LocalEnergyResults &rhs);
    //    //static bool decode(const Node& node, MolecularGeometry& rhs);
    //};
    Emitter& operator<< (Emitter& out, const LocalParticleEnergiesCalculator::LocalEnergyResults& rhs);


    //template<>
    //struct convert<LocalParticleEnergiesCalculator> {
    //    static Node encode(const LocalParticleEnergiesCalculator &rhs);
    //};
    Emitter& operator<< (Emitter& out, const LocalParticleEnergiesCalculator& rhs);

}


#endif //INPSIGHTS_LOCALPARTICLEENERGIESCALCULATION_H
