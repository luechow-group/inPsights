//
// Created by Moria on 19.05.2017.
//

#include "TestElectronAssigner.h"

Assignation TestElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Particle> &electrons) {
    int currentElectron=0;
    Assignation toReturn;
    for(int i=0;i<cores.size();i++){
        toReturn.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            toReturn[i].second.push_back(currentElectron++);
        }
    }
    return toReturn;
}

Assignation TestElectronAssigner::assign(const std::vector<Core> &cores, const std::vector<Electron> &electrons) {
    int currentElectron=0;
    Assignation toReturn;
    for(int i=0;i<cores.size();i++){
        toReturn.emplace_back(i,std::vector<int>());
        for(int j=0;j<cores[i].getCharge();j++){
            toReturn[i].second.push_back(currentElectron++);
        }
    }
    return toReturn;
}
