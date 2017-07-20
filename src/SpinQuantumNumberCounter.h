//
// Created by Morian Sonneton 23.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
#define LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
#include<vector>
#include<utility>


class SpinQuantumNumberCounter {
public:
    SpinQuantumNumberCounter();
    virtual ~SpinQuantumNumberCounter();
    const int & addNumber(const int &newNumber);
    const std::vector<std::tuple<int,int,int> >& getSpinQuantumNumbers() const;
    int getCountOfSpinQuantumNumber(int SpinQuantumNumber) const;
    void printStatsSpinQuantumNumber();
    void printStatsMultiplicities();
    const std::vector<std::tuple<int,int,double> > &getMultiplicities();
private:
    void sortSpinQuantumNumbers();
    static bool pairCompare(const std::tuple<int,int,int> &pair1, const std::tuple<int,int,int> &pair2);

    /*
     *  SpinQuantumNumbers datastructure
     *  Each tuple represents <abs(Ms),count of positive Ms, count of negative Ms>
     *  for Ms=0 the tuple represents <0, count of Ms=0, 0>
     */
    std::vector<std::tuple<int,int,int> > SpinQuantumNumbers;

    /*
     *  Multiplicities datastructure
     *  Each tuple represents <2*S,count of S,relative count of S in percent>
     */
    std::vector<std::tuple<int,int,double> > Multiplicities;
};


#endif //LOCALSPINMULTIPLICITY_SPINQUANTUMNUMBERCOUNTER_H
