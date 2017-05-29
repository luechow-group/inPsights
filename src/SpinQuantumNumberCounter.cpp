//
// Created by Moria on 23.05.2017.
//

#include "SpinQuantumNumberCounter.h"
#include <algorithm>
#include <iostream>
#include <tuple>

bool SpinQuantumNumberCounter::pairCompare(const std::tuple<int, int, int> &tuple1, const std::tuple<int, int, int> &tuple2) {
    return std::get<0>(tuple1)>std::get<0>(tuple2);
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

/*
 * function sorting SpinQuantumNumbers: big abs(Ms) to small abs(Ms)
 */
void SpinQuantumNumberCounter::sortSpinQuantumNumbers() {
    std::sort(SpinQuantumNumbers.begin(),SpinQuantumNumbers.end(),SpinQuantumNumberCounter::pairCompare);
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
    this->sortSpinQuantumNumbers();
    std::cout << "Statistics Spin Quantum Number\n"
            "SpinQuantumNumber Count" << std::endl;
    for(std::vector<std::tuple<int,int,int> >::iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){
        std::cout << static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<1>(*i) << std::endl;
        if(std::get<0>(*i)!=0) std::cout << -static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<2>(*i) << std::endl;
    }
}

/*
 * The count of Mutplicity S is given by adding the values for Ms=-S to +S.
 * In theory for a given S value all Ms values -S to S would occur in same amount, leading to a symmetric distribution over the Ms values
 * Practically, one can see a skewed distribution which is due to interpretation of Multiplicity as Ms value in amolqc.
 * This is why the counts of a particular Ms value belonging to a specific S value is evaluated using linear interpolation.
 */
const std::vector<std::tuple<int,int,double> > &SpinQuantumNumberCounter::getMultiplicities() {
    if(!SpinQuantumNumbers.size())return Multiplicities;
    this->sortSpinQuantumNumbers();                     //sorting of SpinQuantumNumbers is crucial for following algorithm
    for(std::vector<std::tuple<int,int,int> >::iterator i=SpinQuantumNumbers.begin();i!=SpinQuantumNumbers.end();i++){          //i is iterating over the S values
        Multiplicities.emplace_back(std::get<0>(*i),std::get<1>(*i)+std::get<2>(*i),0);
        const int &Msi=std::get<0>(*i);
        const int &countiPositive=std::get<1>(*i);
        const int &countiNegative=std::get<2>(*i);
        for(std::vector<std::tuple<int,int,int> >::iterator j=i+1;j!=SpinQuantumNumbers.end();j++){                             //j is iterationg over the Ms values with abs(Ms)<S
            const int &Msj=std::get<0>(*j);
            int &countjPositive=std::get<1>(*j);
            int &countjNegative=std::get<2>(*j);
            int pSubCounts=((Msi-Msj)*countiNegative+(Msi+Msj)*countiPositive)/(2*Msi);
            int nSubCounts=((Msi-Msj)*countiPositive+(Msi+Msj)*countiNegative)/(2*Msi);
            pSubCounts=pSubCounts>countjPositive?countjPositive:pSubCounts;                                                     //sometimes the count of particular Ms value is smaller than the expected value by linear interpolation.
            nSubCounts=nSubCounts>countjNegative?countjNegative:nSubCounts;                                                     //In such situation only the real counts are taken into account
            std::get<1>(*Multiplicities.rbegin())+=pSubCounts+nSubCounts;
            countjPositive-=pSubCounts;
            countjNegative-=nSubCounts;
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

