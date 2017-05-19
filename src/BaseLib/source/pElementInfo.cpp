#include "pElementInfo.h"

namespace Elements {

    namespace internal {

        pElementInfo* pElementInfo::d_instance = nullptr;

        pElementInfo::pElementInfo() {
            init_data();
        }

        const pElementInfo::ElementData& pElementInfo::operator[] (const ElementType& type) const {
            return d_container.at(type);
        }

        const pElementInfo::ValueType& pElementInfo::operator[] (const std::string& symbol) const {
            for (auto it = d_container.cbegin(); it != d_container.cend(); ++it) {
                if (it->second.symbol() == symbol)
                    return *it;
            }
            throw ElementSymbolNotFound();
        }

        const pElementInfo& pElementInfo::instance() {
            if (d_instance == 0)
                d_instance = new pElementInfo();
            return *d_instance;
        }

        void pElementInfo::init_data() {
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::H,  ElementData("H" ,   1, 1.0079, {237,255,255}, 120.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::He, ElementData("He",   2, 4.0026, {213,255,255}, 140.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Li, ElementData("Li",   3, 6.941 , {204,139,254}, 182.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Be, ElementData("Be",   4, 9.0122, {196,246, 11}, 153.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::B,  ElementData("B" ,   5, 10.811, {255,181,181}, 192.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::C,  ElementData("C" ,   6, 12.011, {144,144,144}, 170.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::N,  ElementData("N" ,   7, 14.007, { 74,112,227}, 155.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::O,  ElementData("O" ,   8, 15.999, {204, 51, 49}, 152.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::F,  ElementData("F" ,   9, 18.988, {148,218,104}, 147.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ne, ElementData("Ne",  10, 20.180, {173,237,244}, 154.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Na, ElementData("Na",  11, 22.990, {168,126,215}, 227.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mg, ElementData("Mg",  12, 24.305, {160,217, 20}, 173.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Al, ElementData("Al",  13, 26.982, {227,161,160}, 184.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Si, ElementData("Si",  14, 28.086, {240,200,160}, 210.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::P,  ElementData("P" ,  15, 30.974, {255,128, 00}, 180.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::S,  ElementData("S" ,  16, 32.065, {231,247, 34}, 180.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cl, ElementData("Cl",  17, 35.453, {105,238, 42}, 175.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ar, ElementData("Ar",  18, 39.948, {139,215,227}, 188.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::K,  ElementData("K" ,  19, 39.098, {136,107,180}, 275.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ca, ElementData("Ca",  20, 40.078, {122,190, 24}, 231.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sc, ElementData("Sc",  21, 44.956, {230,230,230} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ti, ElementData("Ti",  22, 47.867, {191,194,199} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::V,  ElementData("V" ,  23, 50.942, {166,166,171} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cr, ElementData("Cr",  24, 51.996, {138,153,199} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mn, ElementData("Mn",  25, 54.938, {156,122,199} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fe, ElementData("Fe",  26, 55.938, {224,102, 51} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Co, ElementData("Co",  27, 58.933, {240,144,160} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ni, ElementData("Ni",  28, 58.693, { 80,208, 80}, 163.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cu, ElementData("Cu",  29, 63.546, {200,128, 51}, 140.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Zn, ElementData("Zn",  30, 65.38 , {125,128,176}, 139.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ga, ElementData("Ga",  31, 69.723, {204,138,136}, 187.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ge, ElementData("Ge",  32, 72.64 , {154,161,147}, 211.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::As, ElementData("As",  33, 74.922, {189,128,227}, 185.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Se, ElementData("Se",  34, 78.96 , {234,168, 18}, 190.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Br, ElementData("Br",  35, 79.904, {150, 57, 41}, 185.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Kr, ElementData("Kr",  36, 83.798, {109,191,207}, 202.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rb, ElementData("Rb",  37, 83.468, {108, 84,149}, 303.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sr, ElementData("Sr",  38, 87.62 , { 83,165, 24}, 249.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Y,  ElementData("Y" ,  39, 88.906, {135,255,255} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Zr, ElementData("Zr",  40, 91.224, {117,234,234} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Nb, ElementData("Nb",  41, 92.906, { 98,213,215} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mo, ElementData("Mo",  42, 95.96 , { 79,192,196} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tc, ElementData("Tc",  43, 98.91 , { 60,171,179} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ru, ElementData("Ru",  44, 101.07, { 40,150,163} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rh, ElementData("Rh",  45, 102.91, { 20,128,148} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pd, ElementData("Pd",  46, 106.42, {  0,107,134}, 163.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ag, ElementData("Ag",  47, 107.87, {192,192,192}, 172.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cd, ElementData("Cd",  48, 112.41, {255,217,143}, 158.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::In, ElementData("In",  49, 114.82, {186,112,108}, 193.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sn, ElementData("Sn",  50, 118.71, {101,125,126}, 217.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sb, ElementData("Sb",  51, 121.76, {158, 99,181}, 206.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Te, ElementData("Te",  52, 127.60, {208,115, 03}, 206.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::I,  ElementData("I" ,  53, 126.90, {148, 00,148}, 198.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Xe, ElementData("Xe",  54, 131.29, { 81,163,181}, 216.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cs, ElementData("Cs",  55, 132.91, { 85, 56,123}, 343.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ba, ElementData("Ba",  56, 137.33, { 42,142, 20}, 268.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::La, ElementData("La",  57, 138.91, {237,183, 84} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ce, ElementData("Ce",  58, 140.12, {228,187, 83} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pr, ElementData("Pr",  59, 140.91, {221,181, 80} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Nd, ElementData("Nd",  60, 144.24, {214,169, 77} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pm, ElementData("Pm",  61, 146.90, {207,155, 73} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sm, ElementData("Sm",  62, 150.36, {201,140, 68} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Eu, ElementData("Eu",  63, 151.96, {195,126, 64} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Gd, ElementData("Gd",  64, 157.25, {190,112, 59} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tb, ElementData("Tb",  65, 158.93, {184,100, 55} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Dy, ElementData("Dy",  66, 162.50, {179, 89, 51} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ho, ElementData("Ho",  67, 164.93, {173, 79, 48} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Er, ElementData("Er",  68, 167.26, {166, 71, 45} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tm, ElementData("Tm",  69, 168.93, {156, 64, 44} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Yb, ElementData("Yb",  70, 173.05, {142, 60, 45} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Lu, ElementData("Lu",  71, 174.97, {121, 58, 47} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hf, ElementData("Hf",  72, 178.49, {199,183,183} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ta, ElementData("Ta",  73, 180.95, {187,139,174} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::W,  ElementData("W" ,  74, 183.84, {174, 92,162} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Re, ElementData("Re",  75, 186.21, {154, 94,142} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Os, ElementData("Os",  76, 190.23, {133, 97,120} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ir, ElementData("Ir",  77, 192.22, {114, 95,102} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pt, ElementData("Pt",  78, 195.08, {208,208,224}, 175.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Au, ElementData("Au",  79, 196.97, {255,209, 35}, 166.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hg, ElementData("Hg",  80, 200.59, {184,184,208}, 155.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tl, ElementData("Tl",  81, 204.38, {166, 84, 77}, 196.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pb, ElementData("Pb",  82, 207.2 , { 87, 89, 97}, 202.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bi, ElementData("Bi",  83, 208.98, {158, 79,181}, 207.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Po, ElementData("Po",  84, 209.98, {171, 92, 00}, 197.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::At, ElementData("At",  85, 210   , {117, 79, 69}, 202.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rn, ElementData("Rn",  86, 222   , { 56,132,151}, 220.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fr, ElementData("Fr",  87, 223   , { 65, 22,102}, 348.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ra, ElementData("Ra",  88, 226.03, { 00,121, 12}, 283.0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ac, ElementData("Ac",  89, 227   , { 82,183,252} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Th, ElementData("Th",  90, 232.04, { 92,171,240} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pa, ElementData("Pa",  91, 231.04, {101,160,229} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::U,  ElementData("U" ,  92, 238.03, {110,149,218} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Np, ElementData("Np",  93, 237.05, {118,139,208} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pu, ElementData("Pu",  94, 244.10, {126,129,197} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Am, ElementData("Am",  95, 243.10, {133,120,188} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cm, ElementData("Cm",  96, 247.10, {140,111,178} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bk, ElementData("Bk",  97, 247.10, {146,102,169} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cf, ElementData("Cf",  98, 251.10, {152, 94,160} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Es, ElementData("Es",  99, 254.10, {157, 86,151} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fm, ElementData("Fm", 100, 257.10, {162, 78,143} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Md, ElementData("Md", 101, 258   , {166, 71,135} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::No, ElementData("No", 102, 259   , {169, 65,128} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Lr, ElementData("Lr", 103, 262   , {172, 59,120} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rf, ElementData("Rf", 104, 261   , {174, 53,114} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Db, ElementData("Db", 105, 262   , {176, 47,107} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sg, ElementData("Sg", 106, 266   , {178, 42,101} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bh, ElementData("Bh", 107, 264   , {179, 38, 95} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hs, ElementData("Hs", 108, 277   , {179, 34, 90} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mt, ElementData("Mt", 109, 268   , {179, 30, 84} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ds, ElementData("Ds", 110, 281   , {178, 27, 80} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rg, ElementData("Rg", 111, 280   , {177, 24, 75} )));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cn, ElementData("Cn", 112, 285   , {175, 21, 71} )));
        }

    } // namespace internal

} // namespace Elements
