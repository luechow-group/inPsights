//
// Created by Michael Heuer on 2019-02-03.
//

#include <EnergyStatistics.h>


namespace EnergyStatistics {
    ElectronicEnergy::ElectronicEnergy(
            const SingleParticlesStatistics &Te,
            const IntraParticlesStatistics &Vee,
            const InterParticlesStatistics &Ven)
            :
            numberOfElectrons_(Te.mean().size()),
            Te_(Te), Vee_(Vee), Ven_(Ven) {}

    const SingleParticlesStatistics &ElectronicEnergy::Te() const {
        return Te_;
    }

    const IntraParticlesStatistics &ElectronicEnergy::Vee() const {
        return Vee_;
    }

    const InterParticlesStatistics &ElectronicEnergy::Ven() const {
        return Ven_;
    }

    Eigen::Index ElectronicEnergy::numberOfElectrons() const {
        return numberOfElectrons_;
    }
}
