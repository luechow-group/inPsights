#include <utility>

//
// Created by Michael Heuer on 02.10.18.
//

#ifndef INPSIGHTS_ENERGYCALCULATOR_H
#define INPSIGHTS_ENERGYCALCULATOR_H

#include <Statistics.h>
#include <Reference.h>
#include <Sample.h>
#include <Logger.h>
#include <spdlog/spdlog.h>
#include <CoulombPotential.h>
#include <SpinCorrelation.h>

class EnergyCalculator{
public:

    struct TotalEnergies{
        TotalEnergies()
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

    TotalEnergies calculateTotalEnergies() {
        TotalEnergies energies;
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


    unsigned long addReference(const Reference &reference) {
        auto count = unsigned(reference.count());

        Eigen::VectorXd value(1); value[0] = reference.value();
        Eigen::MatrixXd spinCorrelations_ = SpinCorrelation::spinCorrelations(reference.maximum().typesVector()).cast<double>();
        Eigen::VectorXd Te_ = samples_[reference.ownId()].kineticEnergies_;
        Eigen::MatrixXd Vee_ = CoulombPotential::energies(samples_[reference.ownId()].sample_);
        Eigen::MatrixXd Ven_ = CoulombPotential::energies(samples_[reference.ownId()].sample_,atoms_);

        valueStats_.add(value, count);
        spinCorrelationsStats_.add(spinCorrelations_, count);
        TeStats_.add(Te_, count);
        VeeStats_.add(Vee_, count);
        VenStats_.add(Ven_, count);

        return count;
    }

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){
        using namespace YAML;

        yamlDocument_ << BeginDoc << BeginMap
        << Key << "Atoms" << Value << atoms_ << Comment("[a0]")
        << Key << "Vnn" << Comment("[Eh]") << Value << VnnStats_
        << Key << "NSamples" << Value << samples_.size()
        << Key << "Clusters" << Value << BeginSeq;

        size_t totalCount = 0;
        for (auto& cluster : clusteredGloballySimilarMaxima) {

            auto types = cluster[0].representativeReference().maximum().typesVector().asEigenVector();

            for (auto &simRefVector : cluster) {
                valueStats_.reset();
                spinCorrelationsStats_.reset();
                TeStats_.reset();
                VeeStats_.reset();
                VenStats_.reset();

                // Iterate over references being similar to the representative reference.
                std::vector<ElectronsVector> structures;
                for (const auto &ref : simRefVector.similarReferencesIterators()){
                    totalCount += addReference(*ref);
                    structures.push_back(ref->maximum());
                }

                printCluster(structures);
            }
        }
        console->info("overall count {}", totalCount);
        assert(totalCount == samples_.size() && "The total count must match the sample size.");

        yamlDocument_ << EndSeq << EndMap << EndDoc;
        assert(yamlDocument_.good());
    }


    void printCluster(std::vector<ElectronsVector>& structures){
        using namespace YAML;

        yamlDocument_
        << BeginMap
        << Key << "N" << Value << TeStats_.getTotalWeight()
        << Key << "ValueRange" << Value << Comment("[]")
        << valueStats_;

        yamlDocument_ << Key << "Structures" << Comment("[a0]") << Value << BeginSeq;

        size_t nWanted = 16;
        if(structures.size() < nWanted){
            for (auto& i : structures) yamlDocument_ << i;
        } else {
            size_t skipLength = (structures.size()-1) / (nWanted-1);
            for (size_t i = 0; i < nWanted; ++i) {
                yamlDocument_ << structures[i*skipLength];
            }
        }

        yamlDocument_
        << EndSeq << Newline
        << Key << "SpinCorrelations" << Comment("[]") 
        << spinCorrelationsStats_;

        yamlDocument_
        << Key << "Te" << Comment("[Eh]")
        << TeStats_;

        yamlDocument_
        << Key << "Vee" << Comment("[Eh]")
        << VeeStats_;

        yamlDocument_
        << Key << "Ven" << Comment("[Eh]")
        << VenStats_
        << EndMap;
    }

    YAML::Node getYamlNode(){ return YAML::Load(yamlDocument_.c_str()); }

    std::string getYamlDocumentString(){ return std::string(yamlDocument_.c_str()); }

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned,true> spinCorrelationsStats_, VeeStats_, VnnStats_;
    Statistics::RunningStatistics<Eigen::VectorXd,unsigned> TeStats_, valueStats_;
    Statistics::RunningStatistics<Eigen::MatrixXd,unsigned> VenStats_;

    Eigen::MatrixXd Vnn_;

    YAML::Emitter yamlDocument_;
    std::shared_ptr<spdlog::logger> console;
};

#endif //INPSIGHTS_ENERGYCALCULATOR_H
