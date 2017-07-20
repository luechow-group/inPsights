//
// Created by Morian Sonnet on 23.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
#define LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
#include<vector>
#include<utility>

/*
 * This class represents the SpinProjectionQuantumNumberCounter.
 * It is responsible for counting the occurence of different SpinProjectionQuantumNumbers
 * and analyzing the SpinProjectionQuantumNumber-Statistics, returning the occurence of the different
 * SpinProjectionQuantumNumbers.
 */
class SpinProjectionQuantumNumberCounter {
public:
    SpinProjectionQuantumNumberCounter();
    virtual ~SpinProjectionQuantumNumberCounter();
    const int & addNumber(const int &newNumber);
    const std::vector<std::tuple<int,int,int> >& getSpinQuantumNumbers() const;
    int getCountOfSpinProjectionQuantumNumber(int SpinProjectionQuantumNumber) const;
    void printStatsSpinProjectionQuantumNumber();
    void printStatsMultiplicity();
    const std::vector<std::tuple<int,int,double> > &getSpinQuantumNumbers();
private:
    void sortSpinQuantumNumbers();
    static bool pairCompare(const std::tuple<int,int,int> &pair1, const std::tuple<int,int,int> &pair2);

    /*
     *  SpinProjectionQuantumNumbers datastructure
     *  Each tuple represents <abs(Ms),count of positive Ms, count of negative Ms>
     *  for Ms=0 the tuple represents <0, count of Ms=0, 0>
     */
    std::vector<std::tuple<int,int,int> > SpinProjectionQuantumNumbers;

    /*
     *  SpinQuantumNumbers datastructure
     *  Each tuple represents <2*S,count of S,relative count of S in percent>
     */
    std::vector<std::tuple<int,int,double> > SpinQuantumNumbers;
};


#endif //LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
