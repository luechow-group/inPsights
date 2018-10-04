//
// Created by Michael Heuer on 17.08.18.
//

#ifndef AMOLQCPP_COULOMBPOTENTIAL_H
#define AMOLQCPP_COULOMBPOTENTIAL_H

#include <Eigen/Core>
#include "ParticlesVector.h"
#include "Metrics.h"
#include "NaturalConstants.h"
#include "EigenYamlConversion.h"

namespace CoulombPotential {

    template<typename Type>
    Eigen::MatrixXd energies(const ParticlesVector<Type> &pv, bool atomicUnits = true){
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(pv.numberOfEntities(),pv.numberOfEntities());

        for (size_t i = 0; i < V.rows(); i++)
            for (size_t j = i + 1; j < V.cols(); j++)
                V(i,j) = pv[i].charge() * pv[j].charge() / Metrics::distance(pv[i].position(), pv[j].position());

        // symmetrization
        V = V.selfadjointView<Eigen::Upper>();

        // optionally calculate in SI units
        if (!atomicUnits) V * std::pow(Constant::elementaryCharge,2)/(4.0*M_PI*Constant::electricConstant);

        return V;
    };

    template<typename Type1,typename Type2>
    Eigen::MatrixXd energies(const ParticlesVector<Type1> &pv1,
                             const ParticlesVector<Type2> &pv2, bool atomicUnits = true){
        Eigen::MatrixXd V(pv1.numberOfEntities(),pv2.numberOfEntities());

        for (size_t i = 0; i < V.rows(); i++)
            for (size_t j = 0; j < V.cols(); j++)
                V(i,j) = pv1[i].charge() * pv2[j].charge() / Metrics::distance(pv1[i].position(), pv2[j].position());

        // optionally calculate in SI units
        if (!atomicUnits) V * (std::pow(Constant::elementaryCharge,2) /(4.0*Constant::pi*Constant::electricConstant));

        return V;
    };

    template<typename Type>
    void yamlFormattedEnergies(YAML::Emitter & out, const ParticlesVector<Type> &pv) {

        out << YAML::BeginSeq;

        Eigen::MatrixXd V = CoulombPotential::energies(pv);

        for (int i = 0; i < V.rows()-1; ++i)
            out << YAML::Value << Eigen::VectorXd(V.row(i).segment(i+1,V.cols()-(i+1)));

        out << YAML::EndSeq;
    }

    void yamlFormattedEnergies(YAML::Emitter & out, const Eigen::MatrixXd& V, bool noSelfInteractionsQ = false) {
        out << YAML::BeginSeq;

        for (Eigen::Index i = 0; i < V.rows()- (noSelfInteractionsQ? 1 : 0); ++i) {
            if(noSelfInteractionsQ)
                out << YAML::Value << Eigen::VectorXd(V.row(i).segment(i+1,V.cols()-(i+1)));
            else
                out << YAML::Value << Eigen::VectorXd(V.row(i));
        }
        out << YAML::EndSeq;
    }

};

#endif //AMOLQCPP_COULOMBPOTENTIAL_H
