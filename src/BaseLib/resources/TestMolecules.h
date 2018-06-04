//
// Created by Michael Heuer on 30.05.18.
//

#ifndef AMOLQCPP_TESTMOLECULES_H
#define AMOLQCPP_TESTMOLECULES_H

#include <MolecularGeometry.h>

namespace TestMolecules {

    namespace twoElectrons {
        const MolecularGeometry oppositeSpin = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, {0, 0, 0.37}},
                                 {Spin::beta,  {0, 0,-0.37}}})};

        const MolecularGeometry oppositeSpinReversedOrder = {
                AtomsVector(),
                ElectronsVector({{Spin::beta,  {0, 0, 0.37}},
                                 {Spin::alpha, {0, 0,-0.37}}})};

        const MolecularGeometry sameSpinAlpha= {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, {0, 0, 0.37}},
                                 {Spin::alpha, {0, 0,-0.37}}})};
        const MolecularGeometry sameSpinBeta = {
                AtomsVector(),
                ElectronsVector({{Spin::beta, {0, 0, 0.37}},
                                 {Spin::beta, {0, 0,-0.37}}})};
    }
    namespace threeElectrons {
        const MolecularGeometry normal = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, {0, 0, 0.37}},
                                 {Spin::alpha, {0, 0, 0.37}},
                                 {Spin::beta,  {0, 0,-0.37}}})};
    }

    namespace H2 {
        namespace ElectronsInCores {
            const MolecularGeometry normal = {
                    AtomsVector({{Element::H, {0, 0, 0.37}},
                                 {Element::H, {0, 0, -0.37}}}),
                    ElectronsVector({{Spin::alpha, {0, 0, 0.37}},
                                     {Spin::beta,  {0, 0,-0.37}}})};

            const MolecularGeometry translated = {
                    AtomsVector({{Element::H, {0+1, 0+1, 0.37+1}},
                                 {Element::H, {0+1, 0+1, -0.37+1}}}),
                    ElectronsVector({{Spin::alpha, {0+1, 0+1, 0.37+1}},
                                     {Spin::beta,  {0+1, 0+1,-0.37+1}}})};

            const MolecularGeometry reversedElectronOrder = {
                    AtomsVector({{Element::H, {0, 0, 0.37}},
                                 {Element::H, {0, 0,-0.37}}}),
                    ElectronsVector({{Spin::beta,  {0, 0, 0.37}},
                                     {Spin::alpha, {0, 0,-0.37}}})};

            const MolecularGeometry flippedSpins = {
                    AtomsVector({{Element::H, {0, 0, 0.37}},
                                 {Element::H, {0, 0, -0.37}}}),
                    ElectronsVector({{Spin::alpha, {0, 0, -0.37}},
                                     {Spin::beta,  {0, 0, 0.37}}})};
        }

        namespace ElectronsOutsideCores{
            const MolecularGeometry normal = {
                    AtomsVector({{Element::H, {0, 0, 0.37}},
                                 {Element::H, {0, 0, -0.37}}}),
                    ElectronsVector({{Spin::alpha, {0, 0, 0.2}},
                                     {Spin::beta,  {0, 0,-0.2}}})};

            const MolecularGeometry offCenter = {
                    AtomsVector({{Element::H, {0, 0,  0.37}},
                                 {Element::H, {0, 0, -0.37}}}),
                    ElectronsVector({{Spin::alpha, {0, 0.1, 0.2}},
                                     {Spin::beta,  {0,-0.1,-0.2}}})};

            const MolecularGeometry offCenterRotated90 = {
                    AtomsVector({{Element::H, {0, 0,  0.37}},
                                 {Element::H, {0, 0, -0.37}}}),
                    ElectronsVector({{Spin::alpha, { 0.1, 0, 0.2}},
                                     {Spin::beta,  {-0.1, 0,-0.2}}})};
        }
    }

    namespace HeH {
        namespace ElectronsInCores {
            const MolecularGeometry normal = {
                    AtomsVector({{Element::He,{0, 0, 0.37}},
                                 {Element::H, {0, 0,-0.37}}}),
                    ElectronsVector({{Spin::alpha,{0, 0, 0.37}},
                                     {Spin::alpha,{0, 0, 0.37}},
                                     {Spin::beta, {0, 0,-0.37}}})};
        }
    }

    namespace CO2{
        const MolecularGeometry nuclei = {
                AtomsVector({{Element::C,{0,0, 0}},
                             {Element::O,{0,0, 1}},
                             {Element::O,{0,0,-1}}}),
                ElectronsVector()
        };
    }
}

#endif //AMOLQCPP_TESTMOLECULES_H
