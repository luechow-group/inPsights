#ifndef AMOLQCGUI_ELEMENTINFO_H
#define AMOLQCGUI_ELEMENTINFO_H

#include "ElementType.h"
#include "pElementInfo.h"
#include "ElementColor.h"

namespace Elements {


/*
 * Class ElementInfo:
 * Provides information about elements, such as mass, van-der-Waals radius, etc.
 */
    class ElementInfo {
    public:
        static ElementType elementTypeForSymbol(std::string symbol);

        static std::string symbol(ElementType e);

        static double mass(ElementType e);

        /*! Returns the van der Waals radius in atomic units. */
        static double vdwRadius(ElementType e);

        static unsigned Z(ElementType e);

        /* Number of valence electrons */
        static int valElectrons(ElementType e);
        /* Number of s-valence electrons */
        static int sElectrons(ElementType e);
        /* Number of p-valence electrons */
        static int pElectrons(ElementType e);
        /* Number of d-valence electrons */
        static int dElectrons(ElementType e);
        /* Number of f-valence electrons */
        static int fElectrons(ElementType e);

        /* Element ElememtColor for plotting */
        static ElememtColor color(ElementType e);
    };

} // namespace Elements

#endif // AMOLQCGUI_ELEMENTINFO_H
