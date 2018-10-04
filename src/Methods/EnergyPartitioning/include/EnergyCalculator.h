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
        : Te(0),Vee(0),Ven(0),Vnn(0){};

        double totalEnergy() const {
            return Te + Vee + Ven + Vnn;
        };

        double Te, Vee, Ven, Vnn;
    };

    EnergyCalculator(const std::vector<Sample>& samples, const AtomsVector& atoms)
    :
    atoms_(atoms),
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

        TeStats.add(Te_, unsigned(count));
        VeeStats.add(Vee_, unsigned(count));
        VenStats.add(Ven_, unsigned(count));

        return count;
    }

    //TODO eliminate by permuting the samples
    unsigned long addEnergies(const Reference &reference, const Eigen::PermutationMatrix<Eigen::Dynamic>& perm){
        unsigned long count = reference.count();

        auto sampleCopy = samples_[reference.id_].sample_;
        sampleCopy.permute(perm);

        Te_ = perm* (samples_[reference.id_].kineticEnergies_);
        Vee_ = CoulombPotential::energies(sampleCopy);
        Ven_ = CoulombPotential::energies(sampleCopy,atoms_);

        TeStats.add(Te_, unsigned(count));
        VeeStats.add(Vee_, unsigned(count));
        VenStats.add(Ven_, unsigned(count));

        return count;
    }

    void calculateStatistics(const std::vector<std::vector<SimilarReferences>>& clusteredGloballySimilarMaxima){

        yamlDoc_
                << YAML::BeginDoc /*Atoms + VNN*/
                << YAML::BeginMap; /*ADD representative maximum*/
        yamlDoc_
            << YAML::Key << "Atoms" << YAML::Value << atoms_
            << YAML::Key << "Vnn";
        printVnn();

        size_t totalCount = 0;

        for (auto& cluster : clusteredGloballySimilarMaxima) {
            for (auto &simRefVector : cluster) {

                TeStats.reset();
                VeeStats.reset();
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

        yamlDoc_ << YAML::EndMap
                 << YAML::EndDoc;
        std::cout << yamlDoc_.c_str() << std::endl;


    }

    void printTe(){
        yamlDoc_
        << YAML::BeginMap
        << YAML::Key << "mean" << YAML::Newline << YAML::Value << TeStats.mean()
        << YAML::Key << "error"<< YAML::Newline << YAML::Value << TeStats.standardError()
        << YAML::EndMap;
    }

    void printVee(){
        yamlDoc_
        << YAML::BeginMap
        << YAML::Key << "mean" << YAML::Newline
        << YAML::BeginSeq;

        for (int i = 0; i < VeeStats.mean().rows(); ++i) {
            Eigen::VectorXd temp = VeeStats.mean().row(i).segment(i+1,Vee_.cols()-(i+1));
            yamlDoc_ << YAML::Value << temp ;
        }
        yamlDoc_
        << YAML::EndSeq

        << YAML::Key << "error" << YAML::Newline
        << YAML::BeginSeq;
        for (int i = 0; i < VeeStats.standardError().rows(); ++i) {
            Eigen::VectorXd temp = VeeStats.standardError().row(i).segment(i+1,Vee_.cols()-(i+1));
            yamlDoc_ << YAML::Value << temp ;
        }
        yamlDoc_
        << YAML::EndSeq
        << YAML::EndMap;
    }

    void printVen(){
        yamlDoc_
                << YAML::BeginMap
                << YAML::Key << "mean" << YAML::Newline
                << YAML::BeginSeq;

        for (int i = 0; i < VenStats.mean().rows(); ++i) {
            Eigen::VectorXd temp = VenStats.mean().row(i);
            yamlDoc_ << YAML::Value << temp ;
        }
        yamlDoc_
        << YAML::EndSeq
        << YAML::Key << "error" << YAML::Newline
        << YAML::BeginSeq;

        for (int i = 0; i < VenStats.standardError().rows(); ++i) {
            Eigen::VectorXd temp = VenStats.standardError().row(i);
            yamlDoc_ << YAML::Value << temp ;
        }
        yamlDoc_
        << YAML::EndSeq
        << YAML::EndMap;
    }

    void printVnn(){
        yamlDoc_ << YAML::BeginSeq;

        Eigen::MatrixXd Vnn = CoulombPotential::energies(atoms_);

        for (int i = 0; i < Vnn.rows(); ++i) {
            Eigen::VectorXd temp = Vnn.row(i).segment(i+1,Vnn.cols()-(i+1));
            yamlDoc_ << YAML::Value << temp ;
        }
        yamlDoc_ << YAML::EndSeq;
    }

    void printCluster(){

        yamlDoc_
        << YAML::Key << "numberOfSamples" << YAML::Value << TeStats.getTotalWeight()
        << YAML::Key << "Te";
        printTe();

        yamlDoc_
        << YAML::Key << "Vee";
        printVee();

        yamlDoc_
        << YAML::Key << "Ven";
        printVen();
    }

private:
    const std::vector<Sample>& samples_;
    AtomsVector atoms_;

    Statistics::RunningStatistics<Eigen::VectorXd> TeStats;
    Statistics::RunningStatistics<Eigen::MatrixXd> VeeStats, VenStats;


    Eigen::VectorXd Te_;
    Eigen::MatrixXd Vee_, Ven_;

    YAML::Emitter yamlDoc_;
    std::shared_ptr<spdlog::logger> console;
};


#endif //AMOLQCPP_ENERGYCALCULATOR_H
