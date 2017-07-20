//
// Created by Morian Sonneton 19.05.2017.
//

#include "HungarianElectronAssigner.h"
#include <iostream>
#include "Hungarian.h"

Assignment HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {

    resultAssignment.clear();
    DistanceMatrix.resize(electrons.size(),electrons.size());
    generateDistanceMatrix(cores,electrons);
    MatchMatrix.resize(electrons.size(),electrons.size());
    findMatching();
    DistanceMatrix.resize(0,0);
    generateAssignment(cores,electrons);
    MatchMatrix.resize(0,0);
    return resultAssignment;
}

void HungarianElectronAssigner::generateDistanceMatrix(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                DistanceMatrix(rownumber,k)=Particle::distance(cores[i],electrons[k]);
            }
            rownumber++;
        }
    }
   // std::cout << DistanceMatrix << std::endl;
}

void HungarianElectronAssigner::findMatching() {
    Hungarian::findMatching(DistanceMatrix,MatchMatrix,MATCH_MIN);
    //std::cout << "Gematche Matrix: " << std::endl;
    //std::cout << "\n\n" << MatchMatrix << std::endl;
}

void HungarianElectronAssigner::generateAssignment(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        resultAssignment.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                if(MatchMatrix(rownumber,k)>0)resultAssignment[i].second.push_back(k);
            }
            rownumber++;
        }
    }
}

Assignment HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {

    resultAssignment.clear();
    DistanceMatrix.resize(electrons.size(),electrons.size());
    generateDistanceMatrix(cores,electrons);

    //std::cout << "\n\n\n";
    //std::cout << DistanceMatrix << std::endl;

    MatchMatrix.resize(electrons.size(),electrons.size());

    findMatching();
   // std::cout << MatchMatrix << std::endl;

    DistanceMatrix.resize(0,0);
    generateAssignment(cores,electrons);
    MatchMatrix.resize(0,0);

    /*
    std::cout << "Assignment:\n";
    for(std::vector<int>::const_iterator i=resultAssignment[0].second.begin();i!=resultAssignment[0].second.end();i++){
        std::cout << *i << ' ';
    }
    std::cout << std::endl;
    */

    return resultAssignment;
}

void HungarianElectronAssigner::generateDistanceMatrix(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                DistanceMatrix(rownumber,k)=Particle::distance(cores[i],electrons[k]);
            }
            rownumber++;
        }
    }
    //std::cout << DistanceMatrix << std::endl;
}

void HungarianElectronAssigner::generateAssignment(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        resultAssignment.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                if(MatchMatrix(rownumber,k)>0)resultAssignment[i].second.push_back(k);
            }
            rownumber++;
        }
    }
}
