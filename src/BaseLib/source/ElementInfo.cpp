#include "ElementInfo.h"
#include "pElementInfo.h"

namespace Elements {

    ElementType ElementInfo::elementTypeForSymbol(std::string symbol) {
        return internal::pElementInfo::instance()[symbol].first;
    }

    std::string ElementInfo::symbol(ElementType e) {
        return internal::pElementInfo::instance()[e].symbol();
    }

    double ElementInfo::mass(ElementType e) {
        return internal::pElementInfo::instance()[e].mass();
    }

    double ElementInfo::vdwRadius(ElementType e) {
        return internal::pElementInfo::instance()[e].vdWRadius();
    }

    unsigned ElementInfo::Z(ElementType e) {
        // TODO: could actually do the conversion directly from e
        return internal::pElementInfo::instance()[e].Z();
    }

    int ElementInfo::valElectrons(ElementType e) {
        return internal::pElementInfo::instance()[e].valElectrons();
    }

    int ElementInfo::sElectrons(ElementType e) {
        return internal::pElementInfo::instance()[e].sElectrons();
    }

    int ElementInfo::pElectrons(ElementType e) {
        return internal::pElementInfo::instance()[e].pElectrons();
    }

    int ElementInfo::dElectrons(ElementType e) {
        return internal::pElementInfo::instance()[e].dElectrons();
    }

    int ElementInfo::fElectrons(ElementType e) {
        return internal::pElementInfo::instance()[e].fElectrons();
    }

    ElementColor ElementInfo::color(ElementType e) {
        return internal::pElementInfo::instance()[e].color();
    }

} // namespace Elements
