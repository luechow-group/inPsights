// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#include "GeneralStatistics.h"
#include <CoulombPotential.h>

GeneralStatistics::Result GeneralStatistics::calculate(
        Cluster &maxima, std::vector<Sample> &samples, const AtomsVector& atoms) {

    Result result;

    auto Vnnmat = CoulombPotential::energies(atoms);
    result.Vnn_ = 0;
    for (int i = 0; i < atoms.numberOfEntities(); ++i)
        for (int j = i + 1; j < atoms.numberOfEntities(); ++j)
            result.Vnn_ += Vnnmat(i, j);

    for (const auto& sample : samples) {
        auto ESample =  result.Vnn_;

        auto Te = sample.kineticEnergies_.sum();
        ESample += Te;
        result.TeStats_.add(Eigen::Matrix<double,1,1>(sample.kineticEnergies_.sum()));

        double Vee = 0;
        auto Veemat = CoulombPotential::energies(sample.sample_);
        for (long i = 0; i < sample.sample_.numberOfEntities(); ++i) {
            for (long j = i + 1; j < sample.sample_.numberOfEntities(); ++j) {
                Vee += Veemat(i, j);
            }
        }
        ESample += Vee;
        result.VeeStats_.add(Eigen::Matrix<double,1,1>(Vee));

        auto Ven = CoulombPotential::energies(sample.sample_, atoms).sum();
        ESample += Ven;
        result.VenStats_.add(Eigen::Matrix<double,1,1>(Ven));
        result.EStats_.add(Eigen::Matrix<double,1,1>(ESample));
    }

    for (const auto& maximum : maxima) {
        assert(maximum.isLeaf() && "The maxima list used as an input should contain only leaves.");
        result.valueStats_.add(Eigen::Matrix<double,1,1>(maximum.representative()->value()));
    }

    return result;
}

namespace YAML {
    Node convert<GeneralStatistics::Result>::encode(const GeneralStatistics::Result &rhs) {
        Node node;

        node["ValueRange"] = rhs.valueStats_;
        node["Etotal"] = rhs.EStats_;
        node["Te"] = rhs.TeStats_;
        node["Vee"] = rhs.VeeStats_;
        node["Ven"] = rhs.VenStats_;
        node["Vnn"] = rhs.Vnn_;
        return node;
    }

    bool convert<GeneralStatistics::Result>::decode(const Node &node, GeneralStatistics::Result &rhs) {

        rhs = GeneralStatistics::Result(//TODO BETTER NAMES
                node["ValueRange"].as<SingleValueStatistics>(), //TODO rename to -ln|Psi|^2
                node["Etotal"].as<SingleValueStatistics>(),
                node["Te"].as<SingleValueStatistics>(),
                node["Ven"].as<SingleValueStatistics>(),
                node["Vee"].as<SingleValueStatistics>(),
                node["Vnn"].as<double>()
        );

        return true;
    }

    Emitter &operator<<(Emitter &out, const GeneralStatistics::Result &rhs) {
        out << BeginMap
            << Key << "ValueRange" << Value << Comment("[]") << rhs.valueStats_
            << Key << "Etotal" << Comment("[Eh]") << Value << rhs.EStats_
            << Key << "Te" << Comment("[Eh]") << Value << rhs.TeStats_
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.VeeStats_
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.VenStats_
            << Key << "Vnn" << Comment("[Eh]") << Value << rhs.Vnn_
            << EndMap;

        return out;
    }
}