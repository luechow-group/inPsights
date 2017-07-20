//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
#define LOCALSPINMULTIPLICITY_HUNGARIANELECTRONASSIGNER_H
#include "ElectronAssigner.h"
#include "Assignment.h"

class HungarianElectronAssigner: public ElectronAssigner {
public:
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) override;
    virtual Assignment assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) override;
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
