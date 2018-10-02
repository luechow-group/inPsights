//
// Created by Michael Heuer on 02.10.18.
//

#ifndef AMOLQCPP_ENERGYCALCULATOR_H
#define AMOLQCPP_ENERGYCALCULATOR_H

#include <Statistics.h>
#include <Reference.h>
#include <Sample.h>
#include <Logger.h>
#include <spdlog/spdlog.h>
#include <CoulombPotential.h>

//Remove
#include <HungarianHelper.h>

class EnergyCalculator{
public:

    struct Energies{
        Energies()
        : T(0),Vee(0),Ven(0),Vnn(0){};

        double totalEnergy() const {
            return T + Vee + Ven + Vnn;
        };

        double T, Vee, Ven, Vnn;
    };

    EnergyCalculator(const std::vector<Sample>& samples, const AtomsVector& atoms)
    :
    samples_(samples),
    console(spdlog::get(Logger::name)){
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };
    }

    Energies calculateTotalEnergies() {
        Energies energies;
        for (auto& sample : samples_) {
            energies.T += sample.kineticEnergies_.sum();

            auto Veemat = CoulombPotential::energies(sample.sample_);
            for (int i = 0; i < sample.sample_.numberOfEntities(); ++i) {
                for (int j = i + 1; j < sample.sample_.numberOfEntities(); ++j) {
                    energies.Vee += Veemat(i, j);
                }
            }

            auto Venmat = CoulombPotential::energies(sample.sample_, atoms_);
            energies.Ven += Venmat.sum();

            auto Vnnmat = CoulombPotential::energies(atoms_);
            for (int i = 0; i < atoms_.numberOfEntities(); ++i) {
                for (int j = i + 1; j < atoms_.numberOfEntities(); ++j) {
                    energies.Vnn += Vnnmat(i, j);
                }
            }
        }
        energies.T /= samples_.size();
        energies.Vee /= samples_.size();
        energies.Ven /= samples_.size();
        energies.Vnn /= samples_.size();

        return energies;
    }


    unsigned long addEnergies(const Reference &reference) {
        unsigned long count = reference.count();

        ekin = samples_[reference.id_].kineticEnergies_;

        epot = CoulombPotential::energies(samples_[reference.id_].sample_);
        EkinStats.add(ekin, unsigned(count));
        EpotStats.add(epot, unsigned(count));
        return count;
    }

    unsigned long addEnergies(const Reference &reference, const Eigen::PermutationMatrix<Eigen::Dynamic>& perm,const Reference &other){
        unsigned long count = reference.count();

        auto sampleCopy = samples_[reference.id_].sample_;
        sampleCopy.permute(perm);

        ekin = perm* (samples_[reference.id_].kineticEnergies_);
        epot = CoulombPotential::energies(sampleCopy);

        EkinStats.add(ekin, unsigned(count));
        EpotStats.add(epot, unsigned(count));
        return count;
    }

    void calculateStatistics(
            const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){

        size_t totalCount = 0;
        int i=0;

        for (auto& cluster : clusteredGloballySimilarMaxima) {
            for (auto &simRefVector : cluster) {
                i++;

                EkinStats.reset();
                EpotStats.reset();

                // Representative reference
                size_t repRefCount = addEnergies(*simRefVector.repRefIt_);

                size_t simRefCount = 0;
                // Iterate over references being similar to the representative reference.
                for (const auto &simRef : simRefVector.similarReferences_) {
                    simRefCount += addEnergies(*simRef.it_, simRef.perm_,*simRefVector.repRefIt_);
            }

                std::cout << "mean: (" << EkinStats.getTotalWeight() << ")\n" << EkinStats.mean().transpose() << std::endl;
                if (simRefCount >= 2)
                    std::cout << "stdv: (" << EkinStats.getTotalWeight() << ")\n"  << EkinStats.standardDeviation().transpose() << std::endl << std::endl;

                totalCount += repRefCount + simRefCount;
            }
        }
        console->info("overall count {}",totalCount);


    }

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;

    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    std::shared_ptr<spdlog::logger> console;
};

#endif //AMOLQCPP_ENERGYCALCULATOR_H
