//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_SPINDETERMINER_H
#define LOCALSPINMULTIPLICITY_SPINDETERMINER_H
#include "Electron.h"


/*
 * This class represents the SpinDeterminer, which determines the spin of an Electron based on its index.
 * The spin of the first n Electrons is alpha, where n is given by the user.
 */
class SpinDeterminer {
public:
    SpinDeterminer(int highestAlphaIndex);
    virtual ~SpinDeterminer();
    spintype determineSpin(int index) const;
private:
    int highestAlphaIndex;
};


#endif //LOCALSPINMULTIPLICITY_SPINDETERMINER_H
