#include <utility>

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

using namespace YAML;

class EnergyCalculator{
public:

    struct Energies{
        Energies()
        : Te(0),Vee(0),Ven(0),Vnn(0){};

        double totalEnergy() const {
            return Te + Vee + Ven + Vnn;
        };

        double Te, Vee, Ven, Vnn;
    };

    EnergyCalculator(const std::vector<Sample>& samples, AtomsVector atoms)
    :
    atoms_(std::move(atoms)),
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
            energies.Te += sample.kineticEnergies_.sum();

            auto Veemat = CoulombPotential::energies(sample.sample_);
            for (int i = 0; i < sample.sample_.numberOfEntities(); ++i)
                for (int j = i + 1; j < sample.sample_.numberOfEntities(); ++j)
                    energies.Vee += Veemat(i, j);

            auto Venmat = CoulombPotential::energies(sample.sample_, atoms_);
            energies.Ven += Venmat.sum();

        }

        auto Vnnmat = CoulombPotential::energies(atoms_);
        for (int i = 0; i < atoms_.numberOfEntities(); ++i)
            for (int j = i + 1; j < atoms_.numberOfEntities(); ++j)
                energies.Vnn += Vnnmat(i, j);

        energies.Te /= samples_.size();
        energies.Vee /= samples_.size();
        energies.Ven /= samples_.size();

        return energies;
    }


    unsigned long addEnergies(const Reference &reference) {
        unsigned long count = reference.count();

        Te_ = samples_[reference.id_].kineticEnergies_;
        Vee_ = CoulombPotential::energies(samples_[reference.id_].sample_);
        Ven_ = CoulombPotential::energies(samples_[reference.id_].sample_,atoms_);

        TeStats_.add(Te_, unsigned(count));
        VeeStats_.add(Vee_, unsigned(count));
        VenStats.add(Ven_, unsigned(count));

        return count;
    }

    //TODO eliminate by permuting the samples
    unsigned long addEnergies(const Reference &reference, const Eigen::PermutationMatrix<Eigen::Dynamic>& perm){
        unsigned long count = reference.count();

        auto sampleCopy = samples_[reference.id_].sample_;
        sampleCopy.permute(perm);

        Te_ = perm* (samples_[reference.id_].kineticEnergies_);
        Vee_ = CoulombPotential::energies<Spin>(sampleCopy);
        Ven_ = CoulombPotential::energies<Spin,Element>(sampleCopy, atoms_);

        TeStats_.add(Te_, unsigned(count));
        VeeStats_.add(Vee_, unsigned(count));
        VenStats.add(Ven_, unsigned(count));

        return count;
    }

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){

        yamlDocument_ << BeginDoc << BeginMap
        << Key << "Atoms" << Value << atoms_
        << Key << "Vnn" << Value;
        CoulombPotential::yamlFormattedEnergies<Element>(yamlDocument_,atoms_);
        yamlDocument_
        << Key << "NSamples" << Value << samples_.size()
        << Key << "Clusters" << Value << BeginSeq;

        size_t totalCount = 0;
        for (auto& cluster : clusteredGloballySimilarMaxima) {
            for (auto &simRefVector : cluster) {

                TeStats_.reset();
                VeeStats_.reset();
                VenStats.reset();

                // Representative reference
                size_t repRefCount = addEnergies(*simRefVector.repRefIt_);

                size_t simRefCount = 0;

                // Iterate over references being similar to the representative reference.
                for (const auto &simRef : simRefVector.similarReferences_)
                    simRefCount += addEnergies(*simRef.it_, simRef.perm_);


                printCluster();
                totalCount += repRefCount + simRefCount;
            }
        }
        console->info("overall count {}", totalCount);
        assert(totalCount == samples_.size() && "The total count must match the sample size.");

        yamlDocument_ << EndSeq << EndMap << EndDoc;
        assert(yamlDocument_.good());
    }

    void printStats(const Statistics::RunningStatistics<Eigen::VectorXd>& stats) {
        yamlDocument_ << BeginMap
        << Key << "Mean" << Value << stats.mean()
        << YAML::Comment(std::to_string(0) + ":" + std::to_string(stats.mean().size()-1))
        << Key << "Error"<< Value << stats.standardError()
        << YAML::Comment(std::to_string(0) + ":" + std::to_string(stats.mean().size()-1))
        << EndMap;
    }
    void printStats(const Statistics::RunningStatistics<Eigen::MatrixXd>& stats, bool selfInteractionsQ = false){
        yamlDocument_ << BeginMap << Key << "Mean" << Value;
        CoulombPotential::yamlFormattedEnergies(yamlDocument_,stats.mean(), selfInteractionsQ);

        yamlDocument_ << Key << "Error" << Value;
        CoulombPotential::yamlFormattedEnergies(yamlDocument_,stats.standardError(), selfInteractionsQ);
        yamlDocument_ << EndMap;
    }

    void printCluster(){
        yamlDocument_
        << BeginMap
        << Key << "N" << Value << TeStats_.getTotalWeight()
        << Key << "Te";
        printStats(TeStats_);

        yamlDocument_
        << Key << "Vee";
        printStats(VeeStats_,true);

        yamlDocument_
        << Key << "Ven";
        printStats(VenStats);
        yamlDocument_ << EndMap;
    }

    Node getYamlNode(){ return Load(yamlDocument_.c_str()); }

    std::string getYamlDocumentString(){ return std::string(yamlDocument_.c_str()); }

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    Statistics::RunningStatistics<Eigen::VectorXd> TeStats_;
    Statistics::RunningStatistics<Eigen::MatrixXd> VeeStats_, VenStats;


    Eigen::VectorXd Te_;
    Eigen::MatrixXd Vee_, Ven_;

    Emitter yamlDocument_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //AMOLQCPP_ENERGYCALCULATOR_H
