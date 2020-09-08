// Copyright (C) 2020 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ENERGYRESULTSBUNDLE_H
#define INPSIGHTS_ENERGYRESULTSBUNDLE_H

#include <Statistics.h>
#include <ErrorHandling.h>

template<typename Type>
struct EnergyResultsBundle {

    void init(const Type& init) {
        E = init;
        Te = init;
        Vee = init;
        Ven = init;
        Vnn = init;
    }

    Type E, Te, Vee, Ven, Vnn;
};


namespace YAML {
    class Emitter;

    template<typename Type>
    Emitter &operator<<(Emitter &out, const EnergyResultsBundle<Type> &rhs) {
        out << BeginMap
            << Key << "E" << Comment("[Eh]") << Value << rhs.E
            << Key << "Te" << Comment("[Eh]") << Value << rhs.Te
            << Key << "Vee" << Comment("[Eh]") << Value << rhs.Vee
            << Key << "Ven" << Comment("[Eh]") << Value << rhs.Ven
            << Key << "Vnn" << Comment("[Eh]") << Value << rhs.Vnn
            << EndMap;
        return out;
    }
}



#endif //INPSIGHTS_ENERGYRESULTSBUNDLE_H
