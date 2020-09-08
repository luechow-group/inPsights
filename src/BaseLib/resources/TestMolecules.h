// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_TESTMOLECULES_H
#define INPSIGHTS_TESTMOLECULES_H

#include <MolecularGeometry.h>
#include "NaturalConstants.h"

namespace TestMolecules {
    template <typename Type>
    Eigen::Vector3d inbetween(const ParticlesVector<Type>& pv,
            std::pair<Eigen::Index, Eigen::Index> indices, double shift = 0.5){
        assert(indices.first >= 0 && indices.first < pv.numberOfEntities());
        assert(indices.second >= 0 && indices.second < pv.numberOfEntities());
        return shift * (pv.positionsVector()[indices.second]
                       - pv.positionsVector()[indices.first]) + pv.positionsVector()[indices.first];
    }

    namespace twoElectrons {
        const MolecularGeometry oppositeSpin = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::alpha, {0, 0, 0.37}},
                    {Spin::beta,  {0, 0,-0.37}}})};

        const MolecularGeometry oppositeSpinReversedOrder = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::beta,  {0, 0, 0.37}},
                    {Spin::alpha, {0, 0,-0.37}}})};

        const MolecularGeometry sameSpinAlpha= {
                AtomsVector(),
                ElectronsVector({
                    {Spin::alpha, {0, 0, 0.37}},
                    {Spin::alpha, {0, 0,-0.37}}})};
        const MolecularGeometry sameSpinBeta = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::beta, {0, 0, 0.37}},
                    {Spin::beta, {0, 0,-0.37}}})};
    }
    namespace threeElectrons {
        const MolecularGeometry normal = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::alpha, {0, 0, 0.37}},
                    {Spin::alpha, {0, 0, 0.0}},
                    {Spin::beta,  {0, 0,-0.37}}})};
        const MolecularGeometry spinFlipped = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::beta,  {0, 0, 0.37}},
                    {Spin::beta,  {0, 0, 0.0}},
                    {Spin::alpha, {0, 0,-0.37}}})};
        const MolecularGeometry ionic = {
                AtomsVector(),
                ElectronsVector({
                    {Spin::alpha, {0, 0, 0.37}},
                    {Spin::alpha, {0, 0,-0.37}},
                    {Spin::beta,  {0, 0,-0.37}}})};
    }

    namespace fourElectrons {
        const MolecularGeometry normal = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, { 0, 1, 0}},
                                 {Spin::alpha, { 0,-1, 0}},
                                 {Spin::beta,  { 1, 0, 0}},
                                 {Spin::beta,  {-1, 0, 0}}})};
        const MolecularGeometry ionic = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, { 1, 0, 0}},
                                 {Spin::alpha, { 0, 0, 0}},
                                 {Spin::beta,  { 0, 0, 0}},
                                 {Spin::beta,  {-1, 0, 0}}})};
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


        const MolecularGeometry cube = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha, { 1, 1, 1}},
                                 {Spin::beta,  { 1, 1,-1}},
                                 {Spin::alpha, {-1,-1, 1}},
                                 {Spin::beta,  {-1,-1,-1}},
                                 {Spin::alpha, { 1,-1,-1}},
                                 {Spin::beta,  { 1,-1, 1}},
                                 {Spin::alpha, {-1, 1,-1}},
                                 {Spin::beta,  {-1, 1, 1}}
                                })};
    }
    namespace sixteenElectrons {
        const MolecularGeometry twoNestedCubes = {
                AtomsVector(),
                ElectronsVector({{Spin::alpha,{ 1, 1, 1}},
                                 {Spin::beta, { 1, 1,-1}},
                                 {Spin::alpha,{-1,-1, 1}},
                                 {Spin::beta, {-1,-1,-1}},
                                 {Spin::alpha,{ 1,-1,-1}},
                                 {Spin::beta, { 1,-1, 1}},
                                 {Spin::alpha,{-1, 1,-1}},
                                 {Spin::beta, {-1, 1, 1}},

                                 {Spin::beta, { 2, 2, 2}},
                                 {Spin::alpha,{ 2, 2,-2}},
                                 {Spin::beta, {-2,-2, 2}},
                                 {Spin::alpha,{-2,-2,-2}},
                                 {Spin::beta, { 2,-2,-2}},
                                 {Spin::alpha,{ 2,-2, 2}},
                                 {Spin::beta, {-2, 2,-2}},
                                 {Spin::alpha,{-2, 2, 2}}
                })};
    }


    namespace H2 {
        const MolecularGeometry nuclei = {
                AtomsVector({
                    {Element::H, {0, 0, 0.37}},
                    {Element::H, {0, 0, -0.37}}}),{}};
        namespace ElectronsInCores {
            const MolecularGeometry normal = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0, 0.37}},
                        {Spin::beta,  {0, 0,-0.37}}})};

            const MolecularGeometry ionicLeft = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0,-0.37}},
                        {Spin::beta,  {0, 0,-0.37}}})};

            const MolecularGeometry ionicRight = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0, 0.37}},
                        {Spin::beta,  {0, 0, 0.37}}})};

            const MolecularGeometry translated = {
                    AtomsVector({
                        {Element::H, {0+1, 0+1, 0.37+1}},
                        {Element::H, {0+1, 0+1, -0.37+1}}}),
                    ElectronsVector({
                        {Spin::alpha, {0+1, 0+1, 0.37+1}},
                        {Spin::beta,  {0+1, 0+1,-0.37+1}}})};

            const MolecularGeometry reversedElectronOrder = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::beta,  {0, 0, 0.37}},
                        {Spin::alpha, {0, 0,-0.37}}})};

            const MolecularGeometry flippedSpins = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0, -0.37}},
                        {Spin::beta,  {0, 0, 0.37}}})};
        }

        namespace ElectronsOutsideCores{
            const MolecularGeometry normal = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0, 0.2}},
                        {Spin::beta,  {0, 0,-0.2}}})};

            const MolecularGeometry offCenter = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, {0, 0.1, 0.2}},
                        {Spin::beta,  {0,-0.1,-0.2}}})};

            const MolecularGeometry offCenterRotated90 = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, { 0.1, 0, 0.2}},
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

    namespace Li2 {
        const AtomsVector nuclei = {
                AtomsVector({
                                    {Element::Li, {0, 0, 2.673 / 2.0}},
                                    {Element::Li, {0, 0, -2.673 / 2.0}}}
                )};

        const MolecularGeometry normal = {
                nuclei,
                ElectronsVector({
                                        {Spin::alpha, nuclei[0].position()},
                                        {Spin::beta,  nuclei[0].position()},
                                        {Spin::alpha, nuclei[1].position()},
                                        {Spin::beta,  nuclei[1].position()},
                                        {Spin::alpha, inbetween(nuclei, {0, 1}) + Eigen::Vector3d{0, 0.5, 0}},
                                        {Spin::alpha, inbetween(nuclei, {0, 1}) - Eigen::Vector3d{0, 0.5, 0}}
                                })};
    }


    namespace H4 {
        const double a = 1.0;

        namespace ring {
            const MolecularGeometry nuclei = {
                    AtomsVector({
                        {Element::H, {a,  0,  0}},
                        {Element::H, {-a, 0,  0}},
                        {Element::H, {0,  a,  0}},
                        {Element::H, {0,  -a, 0}}}),
                    {}
            };

            const MolecularGeometry fourAlpha = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[3]}
                    })
            };
        }

        namespace linear {
            const MolecularGeometry nuclei = {
                    AtomsVector({
                        {Element::H, {0, 0, -1.5 * a}},
                        {Element::H, {0, 0, -0.5 * a}},
                        {Element::H, {0, 0, +0.5 * a}},
                        {Element::H, {0, 0, +1.5 * a}}}),
                    {}
            };

            const MolecularGeometry ionicA = {
                    nuclei.atoms(), // sorted as in particle kit
                    ElectronsVector({  // ____-a0b2-a1__-__b3
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[1]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[3]}
                    })
            };

            const MolecularGeometry ionicNotinParticleKitSystem = {
                    nuclei.atoms(),
                    ElectronsVector({ // ____-a2b0-a1__-__b3
                        {Spin::beta, nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::beta, nuclei.atoms().positionsVector()[3]}
                    })
            };

            const MolecularGeometry ionicAreflected = {
                    nuclei.atoms(), // sorted for as in particle kit
                    ElectronsVector({ // __b3-a1__-a0b2-____
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]}
                    })
            };

            const MolecularGeometry ionicAreflectedAlphaPermuted = {
                    nuclei.atoms(), // sorted for as in particle kit
                    ElectronsVector({ // __b3-a0__-a1b2-____
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},//permuted w.r.t ionicArefleced
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},//permuted w.r.t ionicArefleced
                        {Spin::beta,  nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]}
                    })
            };

            const MolecularGeometry ionicAreflectedBetaPermuted = {
                    nuclei.atoms(),
                    ElectronsVector({ // __b2-a1__-a0b3-____
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]},//permuted w.r.t ionicArefleced
                        {Spin::beta,  nuclei.atoms().positionsVector()[2]} //permuted w.r.t ionicArefleced
                    })
            };


            const MolecularGeometry ionicAreflectedReorderedNumbering = {
                    nuclei.atoms(), // sorted as in particle kit
                    ElectronsVector({ // a0__-__b2-a1b3-____
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[1]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[2]}
                        // => different to ionicAreflected in chemical mode
                    })
            };
            const MolecularGeometry ionicRealMax1Var1 = { // -lnPsi^2 = -5.585
                    nuclei.atoms(),
                    ElectronsVector({ // a0b3-____-__b2-a1__
                         {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                         {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                         {Spin::beta,  nuclei.atoms().positionsVector()[2]},
                         {Spin::beta,  nuclei.atoms().positionsVector()[0]}
                    })
            };

            const MolecularGeometry ionicRealMax1Var2 = { // -lnPsi^2 = -5.585
                    nuclei.atoms(),
                    ElectronsVector({ // __b2-a1__-___-a0b3
                        {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[3]}
                    })
            };

            const MolecularGeometry ionicRealMax2Var1 = { // -lnPsi^2 = -4.589
                    nuclei.atoms(), // sorted for as in particle kit
                    ElectronsVector({ // a0b3-____-a1b2-____
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]}
                    })
            };

            const MolecularGeometry ionicRealMax2Var2 = { // -lnPsi^2 = -4.589
                    nuclei.atoms(), // sorted for as in particle kit
                    ElectronsVector({ // ____-a0b3-____-a1b2
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[3]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[1]}
                    })
            };
        }
    }

    namespace H6 {
        const double a = 1.0;
        namespace ring {
            const MolecularGeometry nuclei = {
                    AtomsVector({
                        {Element::H, {a,        0,                       0}},
                        {Element::H, {a / 2.0,  std::sqrt(3) * a / 2.0,  0}},
                        {Element::H, {-a / 2.0, std::sqrt(3) * a / 2.0,  0}},
                        {Element::H, {-a,       0,                       0}},
                        {Element::H, {-a / 2.0, -std::sqrt(3) * a / 2.0, 0}},
                        {Element::H, {a / 2.0,  -std::sqrt(3) * a / 2.0, 0}},
                                }),
                    {}
            };
            const MolecularGeometry alphaOnly = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[4]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[5]},
                        })
            };
            const MolecularGeometry alphaShuffle = {
                    nuclei.atoms(),
                    ElectronsVector({
                        alphaOnly.electrons()[3],
                        alphaOnly.electrons()[1],
                        alphaOnly.electrons()[2],
                        alphaOnly.electrons()[0],
                        alphaOnly.electrons()[5],
                        alphaOnly.electrons()[4]
                    })
            };
        }

        namespace linear {
            const MolecularGeometry nuclei = {
                    AtomsVector({
                        {Element::H, {0, 0, -2.5 * a}},
                        {Element::H, {0, 0, -1.5 * a}},
                        {Element::H, {0, 0, -0.5 * a}},
                        {Element::H, {0, 0, +0.5 * a}},
                        {Element::H, {0, 0, +1.5 * a}},
                        {Element::H, {0, 0, +2.5 * a}}
                    }),
                    {}
            };

            const MolecularGeometry ionicA = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[0]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[3]},
                        {Spin::alpha,  nuclei.atoms().positionsVector()[4]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[5]},

                    })
            };

            const MolecularGeometry ionicB = {
                    nuclei.atoms(),
                    ElectronsVector({
                        {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[1]},
                        {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[3]},
                        {Spin::alpha,  nuclei.atoms().positionsVector()[5]},
                        {Spin::beta,  nuclei.atoms().positionsVector()[5]},
                        })
            };
        }
    }

    namespace BH3 {
        const double a = 1.15;
        const MolecularGeometry nuclei = {
                AtomsVector({
                    {Element::B,{ 0, 0, 0}},
                    {Element::H,{ 0, 2.*a, 0}},
                    {Element::H,{-std::sqrt(3)*a,-1.*a, 0}},
                    {Element::H,{ std::sqrt(3)*a,-1.*a, 0}},
                    }),
                {}
        };

        const MolecularGeometry normal = {
                nuclei.atoms(),
                ElectronsVector({{Spin::alpha, nuclei.atoms().positionsVector()[0]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[3]},

                                 {Spin::beta , nuclei.atoms().positionsVector()[0]},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,1},0.25)},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,2},0.25)},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,3},0.25)}
                                })};

        const MolecularGeometry ionic = {
                nuclei.atoms(),
                ElectronsVector({
                                 {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[2]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[3]},

                                 {Spin::beta , nuclei.atoms().positionsVector()[0]},
                                 {Spin::beta , nuclei.atoms().positionsVector()[1]},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,2},0.25)},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,3},0.25)}
                                })};

        const MolecularGeometry ionicMinimal = {
                nuclei.atoms(),
                ElectronsVector({
                    ionic.electrons()[2],
                    ionic.electrons()[3],
                    ionic.electrons()[6],
                    ionic.electrons()[7]})};

        const MolecularGeometry ionicRotated = {
                nuclei.atoms(),
                ElectronsVector({
                                 {Spin::alpha, nuclei.atoms().positionsVector()[0]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[3]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[1]},
                                 {Spin::alpha, nuclei.atoms().positionsVector()[2]},

                                 {Spin::beta , nuclei.atoms().positionsVector()[0]},
                                 {Spin::beta , nuclei.atoms().positionsVector()[3]},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,1},0.25)},
                                 {Spin::beta , inbetween(nuclei.atoms(),{0,2},0.25)}
                                })};


        const MolecularGeometry ionicMinimalRotated = {
                nuclei.atoms(),
                ElectronsVector({
                    ionicRotated.electrons()[2],
                    ionicRotated.electrons()[3],
                    ionicRotated.electrons()[6],
                    ionicRotated.electrons()[7]
                })};

        const MolecularGeometry ionicMinimalRotatedPermuted = {
                nuclei.atoms(),
                ElectronsVector({
                    ionicMinimalRotated.electrons()[0],
                    ionicMinimalRotated.electrons()[2],
                    ionicMinimalRotated.electrons()[1],
                    ionicMinimalRotated.electrons()[3]
                })};
    }

    namespace CO2{
        const MolecularGeometry nuclei = {
                AtomsVector({{Element::C,{0,0, 0}},
                             {Element::O,{0,0, 1}},
                             {Element::O,{0,0,-1}}}),{}};

        const MolecularGeometry nucleiPermuted = {
                AtomsVector({{Element::O,{0,0,-1}},
                             {Element::C,{0,0, 0}},
                             {Element::O,{0,0, 1}}}),{}};

        const MolecularGeometry isolatedNuclei = {
                AtomsVector({{Element::C,{0,0, 0}},
                             {Element::O,{0,0, 10}},
                             {Element::O,{0,0,-10}}}),{}};
    }

    namespace Ethane {
        /* 2      7  6   y x
         *  \      \/    |/
         *   0--*--1     *--z
         *  /\      \
         * 3  4      5
         */

        namespace {
            const double a = 0.7623 * 1.89; // Angstrom to Bohr
            const double b = 1.1544 * 1.89;
            const double c = 0.5053 * 1.89;
            const double d = 0.8752 * 1.89;
        };

        const MolecularGeometry nuclei = {
                AtomsVector({{Element::C,{ 0,   0,-a}},
                             {Element::C,{ 0,   0, a}},
                             {Element::H,{ 0, 2*c,-b}},
                             {Element::H,{-d,  -c,-b}},
                             {Element::H,{ d,  -c,-b}},
                             {Element::H,{ 0,-2*c, b}},
                             {Element::H,{ d,   c, b}},
                             {Element::H,{-d,   c, b}}}
                             ),{}};

        const MolecularGeometry covalent = {
                nuclei.atoms(),
                ElectronsVector({
                    // C Cores
                    {Spin::alpha/* 0*/,nuclei.atoms().positionsVector()[0]},
                    {Spin::beta /* 1*/,nuclei.atoms().positionsVector()[0]},
                    {Spin::alpha/* 2*/,nuclei.atoms().positionsVector()[1]},
                    {Spin::beta /* 3*/,nuclei.atoms().positionsVector()[1]},
                    // H Cores
                    {Spin::alpha/* 4*/,nuclei.atoms().positionsVector()[2]},
                    {Spin::alpha/* 5*/,nuclei.atoms().positionsVector()[3]},
                    {Spin::alpha/* 6*/,nuclei.atoms().positionsVector()[4]},
                    {Spin::beta /* 7*/,nuclei.atoms().positionsVector()[5]},
                    {Spin::beta /* 8*/,nuclei.atoms().positionsVector()[6]},
                    {Spin::beta /* 9*/,nuclei.atoms().positionsVector()[7]},
                    // C-C Bond
                    {Spin::alpha/*10*/,inbetween(nuclei.atoms(),{0,1},0.25)},
                    {Spin::beta /*11*/,inbetween(nuclei.atoms(),{1,0},0.25)},
                    // C-H Bonds
                    {Spin::beta /*12*/,inbetween(nuclei.atoms(),{0,2},0.25)},
                    {Spin::beta /*13*/,inbetween(nuclei.atoms(),{0,3},0.25)},
                    {Spin::beta /*14*/,inbetween(nuclei.atoms(),{0,4},0.25)},
                    {Spin::alpha/*15*/,inbetween(nuclei.atoms(),{1,5},0.25)},
                    {Spin::alpha/*16*/,inbetween(nuclei.atoms(),{1,6},0.25)},
                    {Spin::alpha/*17*/,inbetween(nuclei.atoms(),{1,7},0.25)}
                })};
    }

     namespace trans13Butadiene{
         const MolecularGeometry nuclei = {
                 AtomsVector({
                     {Element::C,{3.4788910586743689, -0.20684940670984894, 0}},
                     {Element::C,{1.1473849284854649, 0.75322588178767302, 0}},
                     {Element::C,{-1.1473849284854649, -0.75322588178767302, 0}},
                     {Element::C,{-3.4788910586743689, 0.20684940670984894, 0}},
                     {Element::H,{5.1322690178617263, 0.99225732208302109, 0}},
                     {Element::H,{3.7999932986536997, -2.2284971665514668, 0}},
                     {Element::H,{0.89890485824722222, 2.7898969632384936, 0}},
                     {Element::H,{-0.89890485824722222, -2.7898969632384936, 0}},
                     {Element::H,{-3.7999932986536997, 2.2284971665514668, 0}},
                     {Element::H,{-5.1322690178617263, -0.99225732208302109, 0}}}
                 ),{}};
         const MolecularGeometry realA = {
                 nuclei.atoms(),
                 ElectronsVector({
                     {Spin::alpha, {3.5975394223425341, -1.0922217028608261, 0}},
                     {Spin::alpha, {1.1473849284854649, 0.75322588178767302, 0}},
                     {Spin::alpha, {-1.1473849284854649, -0.75322588178767302, 0}},
                     {Spin::alpha, {1.5874049087758011, 0.56391334855723474, 0.55747429951044214}},
                     {Spin::alpha, {5.1322690178617263, 0.99225732208302109, 0}},
                     {Spin::alpha, {-5.1322690178617263, -0.99225732208302109, 0}},
                     {Spin::alpha, {-3.4788910586743689, 0.20684940670984894, 0}},
                     {Spin::alpha, {-0.89890485824722222, -2.7898969632384936, 0}},
                     {Spin::alpha, {0.43047536631193267, 0.35523558781869546, 0}},
                     {Spin::alpha, {-3.7999932986536997, 2.2284971665514668, 0}},
                     {Spin::alpha, {-3.0641941392487668, 0.047909989216852508, -0.53609398440956679}},
                     {Spin::alpha, {1.58734652638192, 0.5639136557992902, -0.55750495427352309}},
                     {Spin::alpha, {3.4788910586743689, -0.20684940670984894, 0}},
                     {Spin::alpha, {0.89890485824722222, 2.7898969632384936, 0}},
                     {Spin::alpha, {-3.0640884005936391, 0.047864782655067044, 0.53606890907805838}},
                     {Spin::beta , {3.4788910586743689, -0.20684940670984894, 0}},
                     {Spin::beta , {4.2008444521336754, 0.26522605974264218, 0}},
                     {Spin::beta , {3.0353734462962407, -0.048588148644622402, -0.55181995436501141}},
                     {Spin::beta , {1.1473849284854649, 0.75322588178767302, 0}},
                     {Spin::beta , {-1.5976478602451296, -0.54587213051120087, 0.5546926293336879}},
                     {Spin::beta , {3.7999932986536997, -2.2284971665514668, 0}},
                     {Spin::beta , {-1.0182362149358659, -1.656055052867011, -0.00027354400496658066}},
                     {Spin::beta , {-0.46520810748895652, -0.34906711760582337, -0.00042587260365004043}},
                     {Spin::beta , {-4.1412187961881823, -0.27989436750681446, 0}},
                     {Spin::beta , {-3.4788910586743689, 0.20684940670984894, 0}},
                     {Spin::beta , {3.0353197029413272, -0.048591232848019632, 0.55179238333649949}},
                     {Spin::beta , {-1.5991959843106891, -0.54534111918739525, -0.55387729127249119}},
                     {Spin::beta , {1.0265302804015863, 1.6592335578777477, 0}},
                     {Spin::beta , {-1.1473849284854649, -0.75322588178767302, 0}},
                     {Spin::beta , {-3.7999932986536997, 2.2284971665514668, 0}}
                 })};
         const MolecularGeometry realB = {
                 nuclei.atoms(),
                 ElectronsVector({
                     {Spin::alpha, {-3.6263041097609268, 1.0715317368586184, 0.075549761493657455}},
                     {Spin::alpha, {-1.1473849284854649, -0.75322588178767302, 0}},
                     {Spin::alpha, {1.1473849284854649, 0.75322588178767302, 0}},
                     {Spin::beta , {-1.5735733865127051, -0.5540220661203088, -0.56557897938306823}},
                     {Spin::alpha, {-5.1322690178617263, -0.99225732208302109, 0}},
                     {Spin::alpha, {5.1322690178617263, 0.99225732208302109, 0}},
                     {Spin::alpha, {3.4788910586743689, -0.20684940670984894, 0}},
                     {Spin::alpha, {0.89890485824722222, 2.7898969632384936, 0}},
                     {Spin::beta , {-0.46447638339523423, -0.31077590039123443, 0.044350899543922594}},
                     {Spin::alpha, {3.7999932986536997, -2.2284971665514668, 0}},
                     {Spin::alpha, {3.0704734359017176, -0.049929929385832052, 0.54270283382236717}},
                     {Spin::alpha, {-1.6266385983175098, -0.5590916491495278, 0.53793715033229095}},
                     {Spin::alpha, {-3.4788910586743689, 0.20684940670984894, 0}},
                     {Spin::alpha, {-0.89890485824722222, -2.7898969632384936, 0}},
                     {Spin::beta , {3.1157920901670964, -0.065714289763476236, -0.57004697565633655}},
                     {Spin::beta , {-3.4788910586743689, 0.20684940670984894, 0}},
                     {Spin::beta , {-4.1917904322939181, -0.28582899664826483, -0.054072308616183992}},
                     {Spin::beta , {-3.0210224711898372, 0.03291394823008887, 0.54378477002307901}},
                     {Spin::beta , {-1.1473849284854649, -0.75322588178767302, 0}},
                     {Spin::alpha, {1.6175906270413587, 0.55384669971988265, -0.55389189804762873}},
                     {Spin::beta , {-3.7999932986536997, 2.2284971665514668, 0}},
                     {Spin::beta , {1.0075308082140824, 1.6424350565946748, -0.044651102810099784}},
                     {Spin::alpha, {0.4693451865292696, 0.34281051369885146, 0.039986240884414201}},
                     {Spin::beta , {4.1436340413104338, 0.25000546773472332, 0.07651350799855218}},
                     {Spin::beta , {3.4788910586743689, -0.20684940670984894, 0}},
                     {Spin::alpha, {-3.0618763233845998, 0.052583017208785707, -0.57344637664668541}},
                     {Spin::beta , {1.6660623806474439, 0.52707856858357605, 0.51914429428366726}},
                     {Spin::beta , {-1.0265964309474509, -1.6483280487161414, 0.050313447089605758}},
                     {Spin::beta , {1.1473849284854649, 0.75322588178767302, 0}},
                     {Spin::beta , {3.7999932986536997, -2.2284971665514668, 0}}
                 })};
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
