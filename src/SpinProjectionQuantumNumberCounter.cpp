//
// Created by Morian Sonnet on 23.05.2017.
//

#include "SpinProjectionQuantumNumberCounter.h"
#include <algorithm>
#include <iostream>
#include <tuple>

bool SpinProjectionQuantumNumberCounter::tupleMSCompare(const std::tuple<int, int, int> &tupleMS1,
                                                        const std::tuple<int, int, int> &tupleMS2) {
    return std::get<0>(tupleMS1)>std::get<0>(tupleMS2);
}

SpinProjectionQuantumNumberCounter::SpinProjectionQuantumNumberCounter() {

}

SpinProjectionQuantumNumberCounter::~SpinProjectionQuantumNumberCounter() {
    spinProjectionQuantumNumbers.clear();
}

const int & SpinProjectionQuantumNumberCounter::addNumber(const int &newNumber) {
    std::vector<std::tuple<int,int,int> >::iterator i;
    for(i=spinProjectionQuantumNumbers.begin();i!=spinProjectionQuantumNumbers.end();i++) {
        if (newNumber == std::get<0>(*i)){
            std::get<1>(*i)++;
            break;
        } else if (newNumber == -std::get<0>(*i)){
            std::get<2>(*i)++;
            break;
        }
    }
    if(i==spinProjectionQuantumNumbers.end()) {
        spinProjectionQuantumNumbers.emplace_back(abs(newNumber),newNumber>=0?1:0,newNumber>=0?0:1);
    }
    return newNumber;
}

/*
 * function sorting spinProjectionQuantumNumbers: big abs(Ms) to small abs(Ms)
 */
void SpinProjectionQuantumNumberCounter::sortSpinProjectionQuantumNumbers() {
    std::sort(spinProjectionQuantumNumbers.begin(), spinProjectionQuantumNumbers.end(),
              SpinProjectionQuantumNumberCounter::tupleMSCompare);
}

const std::vector<std::tuple<int, int, int> > &SpinProjectionQuantumNumberCounter::getSpinProjectionQuantumNumbers() const{
    return spinProjectionQuantumNumbers;
}

int SpinProjectionQuantumNumberCounter::getCountOfSpinProjectionQuantumNumber(int SpinProjectionQuantumNumber) const {
    for(std::vector<std::tuple<int,int,int> >::const_iterator i=spinProjectionQuantumNumbers.begin();i!=spinProjectionQuantumNumbers.end();i++){
        if(SpinProjectionQuantumNumber==std::get<0>(*i))return std::get<1>(*i);
        if(SpinProjectionQuantumNumber==-std::get<0>(*i))return std::get<2>(*i);
    }
    return -1;
}


void SpinProjectionQuantumNumberCounter::printStatsSpinProjectionQuantumNumber(){
    this->sortSpinProjectionQuantumNumbers();
    std::cout << "Statistics SpinProjectionQuantumNumber\n"
            "SpinProjectionQuantumNumber Count" << std::endl;
    for(std::vector<std::tuple<int,int,int> >::iterator i=spinProjectionQuantumNumbers.begin();i!=spinProjectionQuantumNumbers.end();i++){
        std::cout << static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<1>(*i) << std::endl;
        if(std::get<0>(*i)!=0) std::cout << -static_cast<double>(std::get<0>(*i))/2 << '\t' << std::get<2>(*i) << std::endl;
    }
}

/*
 *
 * The count of SpinQuantumNumber S is given by adding the values for Ms=-S to +S.
 * In theory for a given S value all Ms values -S to S would occur in same amount, leading to a symmetric distribution over the Ms values
 * Practically, one can see a skewed distribution which is due to interpretation of Multiplicity as Ms value in amolqc.
 * It can also be due to statistical errors.
 * This is why the counts of a particular Ms value belonging to a specific S value is evaluated using linear interpolation.
 */
void SpinProjectionQuantumNumberCounter::calcSpinQuantumNumbers() {
    if(spinProjectionQuantumNumbers.size()==0)return;             //nothing to do
    this->sortSpinProjectionQuantumNumbers();                     //sorting of spinProjectionQuantumNumbers is crucial for following algorithm
    for(std::vector<std::tuple<int,int,int> >::iterator i=spinProjectionQuantumNumbers.begin();i!=spinProjectionQuantumNumbers.end();i++){          //i is iterating over the S values
        spinQuantumNumbers.emplace_back(std::get<0>(*i),std::get<1>(*i)+std::get<2>(*i),0);
        const int &Msi=std::get<0>(*i);
        const int &countiPositive=std::get<1>(*i);
        const int &countiNegative=std::get<2>(*i);
        for(std::vector<std::tuple<int,int,int> >::iterator j=i+1;j!=spinProjectionQuantumNumbers.end();j++){                             //j is iterating over the Ms values with abs(Ms)<S
            const int &Msj=std::get<0>(*j);
            if((Msj-Msi)%2==1)continue;
            int &countjPositive=std::get<1>(*j);
            int &countjNegative=std::get<2>(*j);
            int pSubCounts=((Msi-Msj)*countiNegative+(Msi+Msj)*countiPositive)/(2*Msi);
            int nSubCounts=((Msi-Msj)*countiPositive+(Msi+Msj)*countiNegative)/(2*Msi);
            pSubCounts=pSubCounts>countjPositive?countjPositive:pSubCounts;                                                     //sometimes the count of particular Ms value is smaller than the expected value by linear interpolation.
            nSubCounts=nSubCounts>countjNegative?countjNegative:nSubCounts;                                                     //In such situation only the real counts are taken into account
            std::get<1>(*spinQuantumNumbers.rbegin())+=pSubCounts+nSubCounts;
            countjPositive-=pSubCounts;
            countjNegative-=nSubCounts;
        }
    }
    //calculate percentages
    int countsum=0;
    for(std::vector<std::tuple<int,int,double> >::const_iterator i=spinQuantumNumbers.begin();i!=spinQuantumNumbers.end();i++){
        countsum+=std::get<1>(*i);
    }
    for(std::vector<std::tuple<int,int,double> >::iterator i=spinQuantumNumbers.begin();i!=spinQuantumNumbers.end();i++){
        std::get<2>(*i)=static_cast<double>(std::get<1>(*i))/countsum*100;
    }
    return;
}

void SpinProjectionQuantumNumberCounter::printStatsMultiplicity(){
    if(!spinQuantumNumbers.size())this->calcSpinQuantumNumbers();
    std::cout<< "Statistics Multiplicity\n"
            "Multiplicity Count Percentage" << std::endl;
    for(std::vector<std::tuple<int,int,double> >::iterator i=spinQuantumNumbers.begin();i!=spinQuantumNumbers.end();i++){
        std::cout << std::get<0>(*i)+1 << '\t' << std::get<1>(*i) << '\t' << std::get<2>(*i) << std::endl;
    }
}

