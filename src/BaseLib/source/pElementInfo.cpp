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
            double defaultRadius = 200.0;
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::H,  ElementData("H" ,   1, 1.0079, {175,175,175}, 120.0,  1, 1, 0, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::He, ElementData("He",   2, 4.0026, {213,255,255}, 140.0,  2, 2, 0, 0, 0)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Li, ElementData("Li",   3, 6.941 , {204,139,254}, 182.0,  1, 3, 0, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Be, ElementData("Be",   4, 9.0122, {196,246, 11}, 153.0,  2, 4, 0, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::B,  ElementData("B" ,   5, 10.811, {255,181,181}, 192.0,  3, 4, 1, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::C,  ElementData("C" ,   6, 12.011, { 50, 50, 50}, 170.0,  4, 4, 2, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::N,  ElementData("N" ,   7, 14.007, { 74,112,227}, 155.0,  5, 4, 3, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::O,  ElementData("O" ,   8, 15.999, {204, 51, 49}, 152.0,  6, 4, 4, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::F,  ElementData("F" ,   9, 18.988, {148,218,104}, 147.0,  7, 4, 5, 0, 0)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ne, ElementData("Ne",  10, 20.180, {173,237,244}, 154.0,  8, 4, 6, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Na, ElementData("Na",  11, 22.990, {168,126,215}, 227.0,  1, 5, 6, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mg, ElementData("Mg",  12, 24.305, {160,217, 20}, 173.0,  2, 6, 6, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Al, ElementData("Al",  13, 26.982, {227,161,160}, 184.0,  3, 6, 7, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Si, ElementData("Si",  14, 28.086, {240,200,160}, 210.0,  4, 6, 8, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::P,  ElementData("P" ,  15, 30.974, {255,128,  1}, 180.0,  5, 6, 9, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::S,  ElementData("S" ,  16, 32.065, {162,173, 25}, 180.0,  6, 6,10, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cl, ElementData("Cl",  17, 35.453, {105,238, 42}, 175.0,  7, 6,11, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ar, ElementData("Ar",  18, 39.948, {139,215,227}, 188.0,  8, 6,12, 0, 0)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::K,  ElementData("K" ,  19, 39.098, {136,107,180}, 275.0,  1, 7,12, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ca, ElementData("Ca",  20, 40.078, {122,190, 24}, 231.0,  2, 8,12, 0, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sc, ElementData("Sc",  21, 44.956, {230,230,230}, defaultRadius, 3, 8,12, 1, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ti, ElementData("Ti",  22, 47.867, {191,194,199}, defaultRadius, 4, 8,12, 2, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::V,  ElementData("V" ,  23, 50.942, {166,166,171}, defaultRadius, 5, 8,12, 3, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cr, ElementData("Cr",  24, 51.996, {138,153,199}, defaultRadius, 6, 7,12, 5, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mn, ElementData("Mn",  25, 54.938, {156,122,199}, defaultRadius, 7, 8,12, 5, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fe, ElementData("Fe",  26, 55.938, {224,102, 51}, defaultRadius, 8, 8,12, 6, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Co, ElementData("Co",  27, 58.933, {240,144,160}, defaultRadius, 9, 8,12, 7, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ni, ElementData("Ni",  28, 58.693, { 80,208, 80}, 163.0, 10, 8,12, 8, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cu, ElementData("Cu",  29, 63.546, {200,128, 51}, 140.0, 11, 7,12,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Zn, ElementData("Zn",  30, 65.38 , {125,128,176}, 139.0, 12, 8,12,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ga, ElementData("Ga",  31, 69.723, {204,138,136}, 187.0, 13, 8,13,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ge, ElementData("Ge",  32, 72.64 , {154,161,147}, 211.0, 14, 8,14,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::As, ElementData("As",  33, 74.922, {189,128,227}, 185.0, 15, 8,15,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Se, ElementData("Se",  34, 78.96 , {234,168, 18}, 190.0, 16, 8,16,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Br, ElementData("Br",  35, 79.904, {150, 57, 41}, 185.0, 17, 8,17,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Kr, ElementData("Kr",  36, 83.798, {109,191,207}, 202.0, 18, 8,18,10, 0)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rb, ElementData("Rb",  37, 83.468, {108, 84,149}, 303.0,  1, 9,18,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sr, ElementData("Sr",  38, 87.62 , { 83,165, 24}, 249.0,  2,10,18,10, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Y,  ElementData("Y" ,  39, 88.906, {135,255,255}, defaultRadius, 3,10,18,11, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Zr, ElementData("Zr",  40, 91.224, {117,234,234}, defaultRadius, 4,10,18,12, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Nb, ElementData("Nb",  41, 92.906, { 98,213,215}, defaultRadius, 5, 9,18,14, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mo, ElementData("Mo",  42, 95.96 , { 79,192,196}, defaultRadius, 6, 9,18,15, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tc, ElementData("Tc",  43, 98.91 , { 60,171,179}, defaultRadius, 7,10,18,15, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ru, ElementData("Ru",  44, 101.07, { 40,150,163}, defaultRadius, 8, 9,18,17, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rh, ElementData("Rh",  45, 102.91, { 20,128,148}, defaultRadius, 9, 9,18,18, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pd, ElementData("Pd",  46, 106.42, {  0,107,134}, 163.0, 10, 8,18,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ag, ElementData("Ag",  47, 107.87, {192,192,192}, 172.0, 11, 9,18,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cd, ElementData("Cd",  48, 112.41, {255,217,143}, 158.0, 12,10,18,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::In, ElementData("In",  49, 114.82, {186,112,108}, 193.0, 13,10,19,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sn, ElementData("Sn",  50, 118.71, {101,125,126}, 217.0, 14,10,20,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sb, ElementData("Sb",  51, 121.76, {158, 99,181}, 206.0, 15,10,21,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Te, ElementData("Te",  52, 127.60, {208,115, 03}, 206.0, 16,10,22,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::I,  ElementData("I" ,  53, 126.90, {148, 00,148}, 198.0, 17,10,23,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Xe, ElementData("Xe",  54, 131.29, { 81,163,181}, 216.0, 18,10,24,20, 0)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cs, ElementData("Cs",  55, 132.91, { 85, 56,123}, 343.0,  1,11,24,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ba, ElementData("Ba",  56, 137.33, { 42,142, 20}, 268.0,  2,12,24,20, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::La, ElementData("La",  57, 138.91, {237,183, 84}, defaultRadius, 3,12,24,21, 0)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ce, ElementData("Ce",  58, 140.12, {228,187, 83}, defaultRadius, 4,12,24,21, 1)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pr, ElementData("Pr",  59, 140.91, {221,181, 80}, defaultRadius, 5,12,24,20, 3)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Nd, ElementData("Nd",  60, 144.24, {214,169, 77}, defaultRadius, 6,12,24,20, 4)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pm, ElementData("Pm",  61, 146.90, {207,155, 73}, defaultRadius, 7,12,24,20, 5)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sm, ElementData("Sm",  62, 150.36, {201,140, 68}, defaultRadius, 8,12,24,20, 6)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Eu, ElementData("Eu",  63, 151.96, {195,126, 64}, defaultRadius, 9,12,24,20, 7)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Gd, ElementData("Gd",  64, 157.25, {190,112, 59}, defaultRadius,10,12,24,21, 7)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tb, ElementData("Tb",  65, 158.93, {184,100, 55}, defaultRadius,11,12,24,20, 9)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Dy, ElementData("Dy",  66, 162.50, {179, 89, 51}, defaultRadius,12,12,24,20,10)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ho, ElementData("Ho",  67, 164.93, {173, 79, 48}, defaultRadius,13,12,24,20,11)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Er, ElementData("Er",  68, 167.26, {166, 71, 45}, defaultRadius,14,12,24,20,12)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tm, ElementData("Tm",  69, 168.93, {156, 64, 44}, defaultRadius,15,12,24,20,13)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Yb, ElementData("Yb",  70, 173.05, {142, 60, 45}, defaultRadius,16,12,24,20,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Lu, ElementData("Lu",  71, 174.97, {121, 58, 47}, defaultRadius,17,12,24,21,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hf, ElementData("Hf",  72, 178.49, {199,183,183}, defaultRadius,18,12,24,22,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ta, ElementData("Ta",  73, 180.95, {187,139,174}, defaultRadius,19,12,24,23,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::W,  ElementData("W" ,  74, 183.84, {174, 92,162}, defaultRadius,20,12,24,24,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Re, ElementData("Re",  75, 186.21, {154, 94,142}, defaultRadius,21,12,24,25,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Os, ElementData("Os",  76, 190.23, {133, 97,120}, defaultRadius,22,12,24,26,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ir, ElementData("Ir",  77, 192.22, {114, 95,102}, defaultRadius,23,12,24,27,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pt, ElementData("Pt",  78, 195.08, {208,208,224}, 175.0, 24,11,24,29,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Au, ElementData("Au",  79, 196.97, {255,209, 35}, 166.0, 25,11,24,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hg, ElementData("Hg",  80, 200.59, {184,184,208}, 155.0, 26,12,24,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Tl, ElementData("Tl",  81, 204.38, {166, 84, 77}, 196.0, 27,12,25,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pb, ElementData("Pb",  82, 207.2 , { 87, 89, 97}, 202.0, 28,12,26,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bi, ElementData("Bi",  83, 208.98, {158, 79,181}, 207.0, 29,12,27,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Po, ElementData("Po",  84, 209.98, {171, 92, 00}, 197.0, 30,12,28,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::At, ElementData("At",  85, 210   , {117, 79, 69}, 202.0, 31,12,29,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rn, ElementData("Rn",  86, 222   , { 56,132,151}, 220.0, 32,12,30,30,14)));
            \
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fr, ElementData("Fr",  87, 223   , { 65, 22,102}, 348.0,  1,13,30,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ra, ElementData("Ra",  88, 226.03, { 00,121, 12}, 283.0,  2,14,30,30,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ac, ElementData("Ac",  89, 227   , { 82,183,252}, defaultRadius, 3,14,30,31,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Th, ElementData("Th",  90, 232.04, { 92,171,240}, defaultRadius, 4,14,30,32,14)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pa, ElementData("Pa",  91, 231.04, {101,160,229}, defaultRadius, 5,14,30,31,16)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::U,  ElementData("U" ,  92, 238.03, {110,149,218}, defaultRadius, 6,14,30,31,17)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Np, ElementData("Np",  93, 237.05, {118,139,208}, defaultRadius, 7,14,30,31,18)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Pu, ElementData("Pu",  94, 244.10, {126,129,197}, defaultRadius, 8,14,30,30,20)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Am, ElementData("Am",  95, 243.10, {133,120,188}, defaultRadius, 9,14,30,30,21)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cm, ElementData("Cm",  96, 247.10, {140,111,178}, defaultRadius,10,14,30,30,23)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bk, ElementData("Bk",  97, 247.10, {146,102,169}, defaultRadius,11,14,30,30,23)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cf, ElementData("Cf",  98, 251.10, {152, 94,160}, defaultRadius,12,14,30,30,24)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Es, ElementData("Es",  99, 254.10, {157, 86,151}, defaultRadius,13,14,30,30,25)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fm, ElementData("Fm", 100, 257.10, {162, 78,143}, defaultRadius,14,14,30,30,26)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Md, ElementData("Md", 101, 258   , {166, 71,135}, defaultRadius,15,14,30,30,27)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::No, ElementData("No", 102, 259   , {169, 65,128}, defaultRadius,16,14,30,30,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Lr, ElementData("Lr", 103, 262   , {172, 59,120}, defaultRadius,17,14,30,31,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rf, ElementData("Rf", 104, 261   , {174, 53,114}, defaultRadius,18,14,30,32,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Db, ElementData("Db", 105, 262   , {176, 47,107}, defaultRadius,19,14,30,33,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Sg, ElementData("Sg", 106, 266   , {178, 42,101}, defaultRadius,20,14,30,34,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Bh, ElementData("Bh", 107, 264   , {179, 38, 95}, defaultRadius,21,14,30,35,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Hs, ElementData("Hs", 108, 277   , {179, 34, 90}, defaultRadius,22,14,30,36,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mt, ElementData("Mt", 109, 268   , {179, 30, 84}, defaultRadius,23,14,30,37,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ds, ElementData("Ds", 110, 281   , {178, 27, 80}, defaultRadius,24,14,30,38,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Rg, ElementData("Rg", 111, 280   , {177, 24, 75}, defaultRadius,25,14,30,39,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Cn, ElementData("Cn", 112, 285   , {175, 21, 71}, defaultRadius,26,14,30,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Nh, ElementData("Nh", 113, 287   , {175,175,175}, defaultRadius,27,14,31,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Fl, ElementData("Fl", 114, 289   , {175,175,175}, defaultRadius,28,14,32,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Mc, ElementData("Mc", 115, 288   , {175,175,175}, defaultRadius,29,14,33,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Lv, ElementData("Lv", 116, 289   , {175,175,175}, defaultRadius,30,14,34,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Ts, ElementData("Ts", 117, 293   , {175,175,175}, defaultRadius,31,14,35,40,28)));
            d_container.insert(std::pair<ElementType,ElementData>(ElementType::Og, ElementData("Og", 118, 294   , {175,175,175}, defaultRadius,32,14,36,40,28)));
        }
    }// namespace internal
} // namespace Elements
