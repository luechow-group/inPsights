//
// Created by Michael Heuer on 2019-01-07.
//

#ifndef INPSIGHTS_SURFACEDATAGENERATOR_H
#define INPSIGHTS_SURFACEDATAGENERATOR_H

#include <Eigen/Core>
#include <DualMC.h>
#include <cstdint>
#include "Volume.h"
#include "Conversion.h"
#include "SurfaceData.h"

//TODO template uint16_t un d uint8_t
class SurfaceDataGenerator {
public:

    SurfaceDataGenerator(const Volume<uint16_t> & volume)
            : volume_(volume) {}

    SurfaceData computeSurfaceData(float iso = 0.5){
        dualmc::DualMC<uint16_t> builder;
        builder.build(&volume_.data.front(),
                      volume_.dimX, volume_.dimY, volume_.dimZ,
                      static_cast<uint8_t>(iso * std::numeric_limits<uint16_t>::max()),
                      true, false, dualMcVertices_, dualMcQuads_);

        surfaceData_.vertices = Conversion::convertVertices(dualMcVertices_);
        surfaceData_.triangles = Conversion::quadsToTriangles(dualMcQuads_);

        return surfaceData_;
    }

private:
    Volume<uint16_t> volume_; //TODO template
    std::vector<dualmc::Vertex> dualMcVertices_;
    std::vector<dualmc::Quad> dualMcQuads_;
    SurfaceData surfaceData_;
};

#endif //INPSIGHTS_SURFACEDATAGENERATOR_H
