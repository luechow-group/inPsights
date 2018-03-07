//
// Created by Michael Heuer on 05.02.18.
//

#ifndef AMOLQCPP_ABSTRACTCOLLECTION_H
#define AMOLQCPP_ABSTRACTCOLLECTION_H

#include <iostream>
#include "Particle.h"

class AbstractCollection {
public:
    explicit AbstractCollection(unsigned long numberOfEntities = 0)
            :numberOfEntities_(numberOfEntities){};


    void add(){ numberOfEntities_++; };

    unsigned long numberOfEntities() const{
        return numberOfEntities_;
    }

    friend std::ostream &operator<<(std::ostream &os, const AbstractCollection &ac){
        for (unsigned long i = 0; i < ac.numberOfEntities(); i++) {
            auto decimalPlaces = unsigned(std::log10(i+1)+1);
            os << std::string(ParticleFormat::significantDigits+3-decimalPlaces, ' ')  << i+1 << ParticleFormat::separator;
        }
        std::cout << std::endl;
        return os;
    }

protected:
    unsigned long numberOfEntities_;
    void setNumberOfEntietes(unsigned long numberOfEntities){
        numberOfEntities_ = numberOfEntities;
    }
};

#endif //AMOLQCPP_ABSTRACTCOLLECTION_H
