//
// Created by Morian Sonnet on 19.05.2017.
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

Assignment HungarianElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {
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
}

void HungarianElectronAssigner::findMatching() {
    Hungarian::findMatching(DistanceMatrix,MatchMatrix,MATCH_MIN);
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
