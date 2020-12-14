// Copyright (C) 2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

// obtained from https://github.com/dominikwodniok/dualmc

#include <DualMC.h>

namespace dualmc {
    Vertex::Vertex(
            VertexComponentsType x,
            VertexComponentsType y,
            VertexComponentsType z
    ) : x(x), y(y), z(z) {}

    
    Quad::Quad(
            QuadIndexType i0,
            QuadIndexType i1,
            QuadIndexType i2,
            QuadIndexType i3
    ) : i0(i0), i1(i1), i2(i2), i3(i3) {}
}
