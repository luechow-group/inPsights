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

class EnergyCalculator{
public:

    EnergyCalculator(const std::vector<Sample>& samples)
    :
    samples_(samples),
    console(spdlog::get(Logger::name)){
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };
    }


    unsigned long addEnergies(const Reference &reference) {
        unsigned long count = reference.count();

        ekin = samples_[reference.id_].kineticEnergies_;

        epot = CoulombPotential::energies(samples_[reference.id_].sample_);
        EkinStats.add(ekin, unsigned(count));
        EpotStats.add(epot, unsigned(count));
        return count;
    }

    unsigned long addEnergies(const Reference &reference, const Eigen::PermutationMatrix<Eigen::Dynamic>& perm) {
        unsigned long count = reference.count();

        //TODO Test this
        auto sampleCopy = samples_[reference.id_].sample_;
        sampleCopy.permute(perm);
        ekin = perm * (samples_[reference.id_].kineticEnergies_);
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
                    simRefCount += addEnergies(*simRef.it_, simRef.perm_);
                }

                console->info("rep ref + sim ref count {}", repRefCount + simRefCount);
                std::cout << "mean: (" << EkinStats.getTotalWeight() << ")\n" << EkinStats.mean().transpose() << std::endl;
                if (simRefCount >= 2)
                    std::cout << "stdv:\n" << EkinStats.standardDeviation().transpose() << std::endl << std::endl;

                totalCount += repRefCount + simRefCount;
            }
        }
        console->info("overall count {}",totalCount);


    }

private:
    const std::vector<Sample>& samples_;

    Statistics::RunningStatistics<Eigen::VectorXd> EkinStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> EpotStats;

    Eigen::VectorXd ekin;
    Eigen::MatrixXd epot;
    std::shared_ptr<spdlog::logger> console;
};

#endif //AMOLQCPP_ENERGYCALCULATOR_H
