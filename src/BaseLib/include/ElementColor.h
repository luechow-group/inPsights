//
// Created by heuer on 09.12.16.
//

#ifndef AMOLQCPP_ELEMENTCOLOR_H
#define AMOLQCPP_ELEMENTCOLOR_H

namespace Elements {
    struct ElementColor {
        ElementColor(unsigned int R, unsigned int G, unsigned int B) : R(R), G(G), B(B) {}
        unsigned int R, G, B;
    };
}

#endif //AMOLQCPP_ELEMENTCOLOR_H
