//
// Created by heuer on 09.12.16.
//

#ifndef INPSIGHTS_ELEMENTCOLOR_H
#define INPSIGHTS_ELEMENTCOLOR_H

namespace Elements {
    struct ElementColor {
        ElementColor(unsigned int R, unsigned int G, unsigned int B) : R(R), G(G), B(B) {}
        unsigned int R, G, B;
    };
}

#endif //INPSIGHTS_ELEMENTCOLOR_H
