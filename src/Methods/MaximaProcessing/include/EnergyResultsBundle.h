/* Copyright (C) 2020 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

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
