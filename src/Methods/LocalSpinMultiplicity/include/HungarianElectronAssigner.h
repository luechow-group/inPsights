//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H

#include "Core.h"
#include "Assignment.h"

/*
 * This class represents a ElectronAssigner using the Hungarian Method for assigning Electrons to Cores.
 * It is thus derived of the abstract base class ElectronAssigner.
 */
class HungarianElectronAssigner{
public:
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons);
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons);
private:
    void generateDistanceMatrix(const std::vector<Core> &cores, const std::vector<Particle> &electrons);
    void generateDistanceMatrix(const std::vector<Core> &cores, const std::vector<Electron> &electrons);
    void findMatching();
    void generateAssignment(const std::vector<Core> &cores, const std::vector<Particle> &electrons);
    void generateAssignment(const std::vector<Core> &cores, const std::vector<Electron> &electrons);
    Eigen::MatrixXd DistanceMatrix,MatchMatrix;
    Assignment resultAssignment;
};


#endif //LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
