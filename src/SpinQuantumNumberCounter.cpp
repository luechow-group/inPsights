//
// Created by Moria on 23.05.2017.
//

#include "SpinQuantumNumberCounter.h"
#include <algorithm>
#include <iostream>
#include <tuple>

bool SpinQuantumNumberCounter::pairCompare(const std::tuple<int, int, int> &pair1, const std::tuple<int, int, int> &pair2) {
    return std::get<0>(pair1)>std::get<0>(pair2);
}

SpinQuantumNumberCounter::SpinQuantumNumberCounter() {

}

SpinQuantumNumberCounter::~SpinQuantumNumberCounter() {
    SpinQuantumNumbers.clear();
}

const int & SpinQuantumNumberCounter::addNumber(const int &newNumber) {
    std::vector<std::tuple<int,int,int> >::iterator i;
    for(i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++) {
        if (newNumber == std::get<0>(*i)){
            std::get<1>(*i)++;
            break;
        } else if (newNumber == -std::get<0>(*i)){
            std::get<2>(*i)++;
            break;
        }
    }
    if(i==SpinQuantumNumbers.end()) {
        SpinQuantumNumbers.emplace_back(abs(newNumber),newNumber>=0?1:0,newNumber>=0?0:1);
    }
    return newNumber;
}

const std::vector<std::tuple<int, int,int> > &SpinQuantumNumberCounter::sort() {
    std::sort(SpinQuantumNumbers.begin(),SpinQuantumNumbers.end(),SpinQuantumNumberCounter::pairCompare);
    return SpinQuantumNumbers;
}

const std::vector<std::tuple<int, int, int> > &SpinQuantumNumberCounter::getSpinQuantumNumbers() const{
    return SpinQuantumNumbers;
}

int SpinQuantumNumberCounter::getCountOfSpinQuantumNumber(int SpinQuantumNumber) const {
    for(std::vector<std::tuple<int,int,int> >::const_iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){
        if(SpinQuantumNumber==std::get<0>(*i))return std::get<1>(*i);
        if(SpinQuantumNumber==-std::get<0>(*i))return std::get<2>(*i);
    }
    return -1;
}


void SpinQuantumNumberCounter::printStatsSpinQuantumNumber(){
    this->sort();
    std::cout << "Statistics Spin Quantum Number\n"
            "SpinQuantumNumber Count" << std::endl;
    for(std::vector<std::tuple<int,int,int> >::iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){
        std::cout << static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<1>(*i) << std::endl;
        if(std::get<0>(*i)!=0) std::cout << -static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<2>(*i) << std::endl;
    }
}

const std::vector<std::tuple<int,int,double> > &SpinQuantumNumberCounter::getMultiplicities() {
    if(!SpinQuantumNumbers.size())return Multiplicities;
    this->sort();
    for(std::vector<std::tuple<int,int,int> >::iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){
        Multiplicities.emplace_back(std::get<0>(*i),std::get<1>(*i)+std::get<2>(*i),0);
        int &Msi=std::get<0>(*i);
        int &countpi=std::get<1>(*i);
        int &countni=std::get<2>(*i);
        for(std::vector<std::tuple<int,int,int> >::iterator j=i+1;j!=SpinQuantumNumbers.end();j++){
            int &Msj=std::get<0>(*j);
            int &countpj=std::get<1>(*j);
            int &countnj=std::get<2>(*j);
            int pSubCounts=((Msi-Msj)*countni+(Msi+Msj)*countpi)/(2*Msi);
            int nSubCounts=((Msi-Msj)*countpi+(Msi+Msj)*countni)/(2*Msi);
            pSubCounts=pSubCounts>countpj?countpj:pSubCounts;
            nSubCounts=nSubCounts>countnj?countnj:nSubCounts;
            std::get<1>(*Multiplicities.rbegin())+=pSubCounts+nSubCounts;
            countpj-=pSubCounts;
            countnj-=nSubCounts;
        }
    }
    //calculate percentages
    int countsum=0;
    for(std::vector<std::tuple<int,int,double> >::const_iterator i=Multiplicities.begin();i!=Multiplicities.end();i++){
        countsum+=std::get<1>(*i);
    }
    for(std::vector<std::tuple<int,int,double> >::iterator i=Multiplicities.begin();i!=Multiplicities.end();i++){
        std::get<2>(*i)=static_cast<double>(std::get<1>(*i))/countsum*100;
    }
    return Multiplicities;
}

void SpinQuantumNumberCounter::printStatsMultiplicities(){
    if(!Multiplicities.size())this->getMultiplicities();
    std::cout<< "Statistics Multiplicity\n"
            "Multiplicity Count Percentage" << std::endl;
    for(std::vector<std::tuple<int,int,double> >::iterator i=Multiplicities.begin();i!=Multiplicities.end();i++){
        std::cout << std::get<0>(*i)+1 << '\t' << std::get<1>(*i) << '\t' << std::get<2>(*i) << std::endl;
    }
}

