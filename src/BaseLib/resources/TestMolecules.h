//
// Created by Michael Heuer on 30.05.18.
//

#ifndef INPSIGHTS_TESTMOLECULES_H
#define INPSIGHTS_TESTMOLECULES_H

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
                                 {Spin::alpha, {0, 0, 0.0}},
                                 {Spin::beta,  {0, 0,-0.37}}})};
    }

    namespace eightElectrons {
        const MolecularGeometry square = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, { 1, 0, 0}},
                                 {Spin::alpha, { 1, 1, 0}},
                                 {Spin::alpha, { 0, 1, 0}},
                                 {Spin::alpha, {-1, 1, 0}},
                                 {Spin::beta,  {-1, 0, 0}},
                                 {Spin::beta,  {-1,-1, 0}},
                                 {Spin::beta,  { 0,-1, 0}},
                                 {Spin::beta,  { 1,-1, 0}}})};
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
                ElectronsVector()};
        const MolecularGeometry isolatedNuclei = {
                AtomsVector({{Element::C,{0,0, 0}},
                             {Element::O,{0,0, 10}},
                             {Element::O,{0,0,-10}}}),
                ElectronsVector()};
    }

    namespace CoulombPotentialTest{
        const MolecularGeometry HeH = {
                AtomsVector({{Element::He,{0, 0,-1}},
                             {Element::H, {0, 0, 1}}}),
                ElectronsVector({{Spin::alpha,{0, 0,-1}},
                                 {Spin::alpha,{0, 0, 0}},
                                 {Spin::beta, {0, 0,-1}}})};
    }
}

#endif //INPSIGHTS_TESTMOLECULES_H
