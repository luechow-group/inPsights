// Copyright 2019 heuer
// SPDX-License-Identifier: GPL-3.0-or-later

#include "ColorPalette.h"

namespace ColorPalette {
    const QColor& colorFunction(size_t i) {
        while (i >= palette.size())
            i -= palette.size();

        return palette[i];
    }
}