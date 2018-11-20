//
// Created by leonard on 22.02.18.
//

#ifndef INPSIGHTS_POTENTIALPROBLEM_H
#define INPSIGHTS_POTENTIALPROBLEM_H

#include <problem.h>
#include <ParticlesVector.h>

class PotentialProblem: public cppoptlib::Problem<double,Eigen::Dynamic>{
public:

    PotentialProblem();

    explicit PotentialProblem(const AtomsVector & atomsVector);

    AtomsVector getAtoms() const;

    double value(const Eigen::VectorXd & x) override;

    void gradient(const Eigen::VectorXd &x, Eigen::VectorXd &grad) override;

    void hessian(const Eigen::VectorXd &x, Eigen::MatrixXd &hessian) override;

private:
    AtomsVector nuclei_;
    double vPotNuclei_;

    void calculateVPotNuclei();
};

#endif //INPSIGHTS_POTENTIALPROBLEM_H
