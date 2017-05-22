//
// Created by Moria on 19.05.2017.
//

#include "HungarianElectronAssigner.h"
#include <iostream>
#include "Hungarian.h"

Assignation HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {

    resultAssignation.clear();
    DistanceMatrix.resize(electrons.size(),electrons.size());
    generateDistanceMatrix(cores,electrons);
    MatchMatrix.resize(electrons.size(),electrons.size());
    findMatching();
    DistanceMatrix.resize(0,0);
    generateAssignation(cores,electrons);
    MatchMatrix.resize(0,0);
    return resultAssignation;
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

void HungarianElectronAssigner::generateAssignation(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        resultAssignation.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                if(MatchMatrix(rownumber,k)>0)resultAssignation[i].second.push_back(k);
            }
            rownumber++;
        }
    }
}

Assignation HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {

    resultAssignation.clear();
    DistanceMatrix.resize(electrons.size(),electrons.size());
    generateDistanceMatrix(cores,electrons);
    MatchMatrix.resize(electrons.size(),electrons.size());
    findMatching();
    DistanceMatrix.resize(0,0);
    generateAssignation(cores,electrons);
    MatchMatrix.resize(0,0);
    return resultAssignation;
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

void HungarianElectronAssigner::generateAssignation(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {
    int rownumber=0;
    for(int i=0;i<cores.size();i++){
        resultAssignation.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            for(int k=0;k<electrons.size();k++){
                if(MatchMatrix(rownumber,k)>0)resultAssignation[i].second.push_back(k);
            }
            rownumber++;
        }
    }
}
