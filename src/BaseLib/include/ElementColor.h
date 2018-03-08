//
// Created by heuer on 09.12.16.
//

#ifndef AMOLQCGUI_ELEMENTCOLOR_H
#define AMOLQCGUI_ELEMENTCOLOR_H

namespace Elements {
    struct ElementColor {
        ElementColor(unsigned int R, unsigned int G, unsigned int B) : R(R), G(G), B(B) {}
        unsigned int R, G, B;
    };
}

#endif //AMOLQCGUI_ELEMENTCOLOR_H
