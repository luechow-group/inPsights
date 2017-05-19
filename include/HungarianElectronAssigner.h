//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
#include "ElectronAssigner.h"
#include "Assignation.h"

class HungarianElectronAssigner: public ElectronAssigner {
public:
    virtual Assignation assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) override;
private:
    void generateDistanceMatrix(const std::vector<Core> &cores, const std::vector<Particle> &electrons);
    void findMatching();
    void generateAssignation(const std::vector<Core> &cores, const std::vector<Particle> &electrons);
    Eigen::MatrixXd DistanceMatrix,MatchMatrix;
    Assignation resultAssignation;
};


#endif //LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
