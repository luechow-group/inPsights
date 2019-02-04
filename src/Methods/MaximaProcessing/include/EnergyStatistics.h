//
// Created by Michael Heuer on 2019-02-03.
//

#ifndef INPSIGHTS_ENERGYSTATISTICS_H
#define INPSIGHTS_ENERGYSTATISTICS_H

#include "GeneralStatistics.h"

namespace EnergyStatistics {

    class ElectronicEnergy {
    public:
        ElectronicEnergy(const SingleParticlesStatistics &Te,
                         const IntraParticlesStatistics &Vee,
                         const InterParticlesStatistics &Ven);

        Eigen::Index numberOfElectrons() const;

        const SingleParticlesStatistics &Te() const;

        const IntraParticlesStatistics &Vee() const;

        const InterParticlesStatistics &Ven() const;

    private:
        Eigen::Index numberOfElectrons_;
        SingleParticlesStatistics Te_; // kinetic energy
        IntraParticlesStatistics Vee_; // electron-electron interaction energy
        InterParticlesStatistics Ven_; // electron-core interaction energy
    };


}

#endif //INPSIGHTS_ENERGYSTATISTICS_H
