//
// Created by Moria on 19.05.2017.
//

#include "HungarianElectronAssigner.h"
#include <iostream>
#include "Hungarian.h"

Assignation HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> distanceMatrix, resultMatrix;      //Matrix is stored as Column first by default
    distanceMatrix.resize(electrons.size(),electrons.size());
    resultMatrix.resize(electrons.size(),electrons.size());
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                distanceMatrix(rownumber,k)=Particle::distance(cores[i],electrons[k]);
            }
            rownumber++;
        }
    }
    std::cout << distanceMatrix << std::endl;
    Hungarian::findMatching(distanceMatrix,resultMatrix,MATCH_MIN);
    std::cout << "und jetzt die neue Matrix" << std::endl;
    std::cout << "\n\n" << resultMatrix << std::endl;

    rownumber=0;
    Assignation toReturn;
    for(int i=0;i<cores.size();i++){
        toReturn.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                if(resultMatrix(rownumber,k)>0)toReturn[i].second.push_back(k);
            }
            rownumber++;
        }
    }
    return toReturn;
}
