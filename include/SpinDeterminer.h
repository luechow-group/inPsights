//
// Created by Moria on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_SPINDETERMINER_H
#define LOCALSPINMULTIPLICITY_SPINDETERMINER_H
#include "Electron.h"


class SpinDeterminer {
public:
    SpinDeterminer(int highestAlphaIndex);
    virtual ~SpinDeterminer();
    spintype determineSpin(int index) const;
private:
    int highestAlphaIndex;

};


#endif //LOCALSPINMULTIPLICITY_SPINDETERMINER_H
