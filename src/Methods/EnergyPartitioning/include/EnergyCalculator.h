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
    samples_(samples),
    atoms_(std::move(atoms)),
    Vnn_(CoulombPotential::energies<Element>(atoms_)),
    console(spdlog::get(Logger::name))
    {
        if(!console){
            Logger::initialize();
            console = spdlog::get(Logger::name);
        };

        VnnStats_.add(Vnn_);
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

        Te_ = samples_[reference.ownId()].kineticEnergies_;
        Vee_ = CoulombPotential::energies(samples_[reference.ownId()].sample_);
        Ven_ = CoulombPotential::energies(samples_[reference.ownId()].sample_,atoms_);

        TeStats_.add(Te_, unsigned(count));
        VeeStats_.add(Vee_, unsigned(count));
        VenStats_.add(Ven_, unsigned(count));

        return count;
    }

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){
        using namespace YAML;

        yamlDocument_ << BeginDoc << BeginMap
        << Key << "Atoms" << Value << atoms_ << Comment("[a0]")
        << Key << "Vnn" << Comment("[Eh]") << Value;
        VnnStats_.toYaml(yamlDocument_, true);

        yamlDocument_
        << Key << "NSamples" << Value << samples_.size()
        << Key << "Clusters" << Value << BeginSeq;

        size_t totalCount = 0;
        for (auto& cluster : clusteredGloballySimilarMaxima) {
            for (auto &simRefVector : cluster) {
                //CAREFUL: this works only for ethane right now, where clusters have size 1 only
                //TODO permute all cluster elements to match each other (difficult, because order has to be kept)
                TeStats_.reset();
                VeeStats_.reset();
                VenStats_.reset();

                // Iterate over references being similar to the representative reference.
                for (const auto &ref : simRefVector.similarReferences_)
                    totalCount += addEnergies(*ref);

                printCluster();
            }
        }
        console->info("overall count {}", totalCount);
        assert(totalCount == samples_.size() && "The total count must match the sample size.");

        yamlDocument_ << EndSeq << EndMap << EndDoc;
        assert(yamlDocument_.good());
    }


    void printCluster(){
        using namespace YAML;

        yamlDocument_
        << BeginMap
        << Key << "N" << Value << TeStats_.getTotalWeight()
        << Key << "Te" << Comment("[Eh]");
        TeStats_.toYaml(yamlDocument_);

        yamlDocument_
        << Key << "Vee" << Comment("[Eh]");
        VeeStats_.toYaml(yamlDocument_,true);

        yamlDocument_
        << Key << "Ven" << Comment("[Eh]");
        VenStats_.toYaml(yamlDocument_);
        yamlDocument_ << EndMap;
    }

    YAML::Node getYamlNode(){ return YAML::Load(yamlDocument_.c_str()); }

    std::string getYamlDocumentString(){ return std::string(yamlDocument_.c_str()); }

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    Statistics::RunningStatistics<Eigen::VectorXd> TeStats_;
    Statistics::RunningStatistics<Eigen::MatrixXd> VeeStats_, VenStats_, VnnStats_;
    Eigen::VectorXd Te_;
    Eigen::MatrixXd Vee_, Ven_, Vnn_;

    YAML::Emitter yamlDocument_;
    std::shared_ptr<spdlog::logger> console;
};

#endif //AMOLQCPP_ENERGYCALCULATOR_H
