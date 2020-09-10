#ifndef INPSIGHTS_ELEMENTTYPE_H
#define INPSIGHTS_ELEMENTTYPE_H

#include <utility>
#include <map>
#include <ostream>

namespace Elements {
    enum class ElementType {
        none = 0,
        H,                                                                                                                          He,
        Li, Be,                                                                                                 B,  C,  N,  O,  F,  Ne,
        Na, Mg,                                                                                                 Al, Si, P,  S,  Cl, Ar,
        K,  Ca, Sc,                                                         Ti, V,  Cr, Mn, Fe, Co, Ni, Cu, Zn, Ga, Ge, As, Se, Br, Kr,
        Rb, Sr, Y,                                                          Zr, Nb, Mo, Tc, Ru, Rh, Pd, Ag, Cd, In, Sn, Sb, Te, I,  Xe,
        Cs, Ba, La, Ce, Pr, Nd, Pm, Sm, Eu, Gd, Tb, Dy, Ho, Er, Tm, Yb, Lu, Hf, Ta, W,  Re, Os, Ir, Pt, Au, Hg, Tl, Pb, Bi, Po, At, Rn,
        Fr, Ra, Ac, Th, Pa, U,  Np, Pu, Am, Cm, Bk, Cf, Es, Fm, Md, No, Lr, Rf, Db, Sg, Bh, Hs, Mt, Ds, Rg, Cn, Nh, Fl, Mc, Lv, Ts, Og
    };

    ElementType first();

    ElementType last();

    ElementType elementFromInt(int type);

    int elementToInt(ElementType element);
    
} // namespace Elements

using Element = Elements::ElementType;

std::ostream& operator<< (std::ostream& os, const Elements::ElementType& e);

namespace YAML {
    class Node; class Emitter;
    template <typename Type> struct convert;

    template<> struct convert<Element> {
        static Node encode(const Element &rhs);
        static bool decode(const Node &node, Element &rhs);
    };
    Emitter &operator<<(Emitter &out, const Element &e);
}

#endif // INPSIGHTS_ELEMENTTYPE_H
