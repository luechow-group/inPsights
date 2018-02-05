#ifndef AMOLQCGUI_ELEMENTTYPES_H
#define AMOLQCGUI_ELEMENTTYPES_H

//#include <iostream>
#include <utility>
#include <map>

namespace Elements {

enum class ElementType {
  none = 0,
  H,                                                                                                                          He,
  Li, Be,                                                                                                 B,  C,  N,  O,  F,  Ne,
  Na, Mg,                                                                                                 Al, Si, P,  S,  Cl, Ar,
  K,  Ca, Sc,                                                         Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
  Rb, Sr, Y,                                                          Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
  Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
  Fr, Ra, Ac, Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn
};

} // namespace Elements

std::ostream& operator<< (std::ostream& os, const Elements::ElementType& e);

#endif // AMOLQCGUI_ELEMENTTYPES_H
