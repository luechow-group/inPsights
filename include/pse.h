//
// Created by Moria on 16.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_PSE_H
#define LOCALSPINMULTIPLICITY_PSE_H
#include <string>

class Pse{
private:
    static const std::string pse[];
public:
    static int findElement(std::string& elementType);
};


#endif //LOCALSPINMULTIPLICITY_PSE_H
