// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_COLORPALETTE_H
#define INPSIGHTS_COLORPALETTE_H

#include <QColor>
#include <vector>

namespace ColorPalette {
    const std::vector<QColor> palette = {
            {237, 75,  4},
            {254, 165, 59},
            {234, 214, 36},
            {134, 1,   238},
            {251, 32,  118},
            {167, 232, 49},
            {99,  129, 35},
            {70,  243, 62},
            {25,  32,  177},
            {169, 92,  154},
            {231, 47,  194},
            {250, 156, 195},
            {75,  25,  93},
            {54,  123, 187},
            {27,  59,  105},
            {101, 230, 249},
            {24,  68,  27},
            {79,  223, 150},
            {22,  146, 148},
            {157, 187, 230},
            {105, 117, 254},
            {224, 250, 182},
            {18,  152, 45},
            {245, 205, 175},
            {83,  2,   8},
            {145, 76,  15},
            {40,  24,  0}
    };

    const QColor& colorFunction(size_t i);
}

#endif //INPSIGHTS_COLORPALETTE_H
