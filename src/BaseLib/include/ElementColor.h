// Copyright (C) 2016-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ELEMENTCOLOR_H
#define INPSIGHTS_ELEMENTCOLOR_H

namespace Elements {
    struct ElementColor {
        ElementColor(unsigned int R, unsigned int G, unsigned int B) : R(R), G(G), B(B) {}
        unsigned int R, G, B;
    };
}

#endif //INPSIGHTS_ELEMENTCOLOR_H
