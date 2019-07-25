#ifndef INPSIGHTS_PELEMENTINFO_H
#define INPSIGHTS_PELEMENTINFO_H

#include "NaturalConstants.h"
#include "ElementType.h"
#include "ElementColor.h"
#include <string>

namespace Elements {


    namespace internal {

/* Class pElementInfo:
   Provides a mapping of ElementType to ElementData and
   accessing Elements by a symbol string.
   Fast lookup for subscript operator [ElementType]. Throws std::out_of_range.
   Slow lookup for subscript operator [std::string]. Throws pElementInfo::DataNotAvailable.
*/
        class pElementInfo {
        public:
            /* Data type for each element of the periodic table
               Data includes: Mass (a.u.), Z,
            */
            class ElementData {
            public:
                class DataNotAvailable {};

                /* Constructor for default (empty) element */
                ElementData() : d_symbol(""), d_Z(0), d_mass(-1),
                                d_vdWRadius(-1), d_valenceElectrons(-1),
                                d_sElectrons(-1), d_pElectrons(-1), d_dElectrons(-1), d_fElectrons(-1),
                                d_color({0,0,0}) {};

                /* Constructor for adding a new element */
                ElementData(std::string symbol, unsigned int Z, double mass,
                            ElementColor color = {0,0,0},
                            double vdWRadiusInPicometers = -1, int valElectrons = -1, int sElectrons = -1,
                            int pElectrons = -1, int dElectrons = -1, int fElectrons = -1) :
                        d_symbol(symbol), d_Z(Z), d_mass(mass),
                        d_vdWRadius(vdWRadiusInPicometers / 100 * ConversionFactors::angstrom2bohr),
                        d_valenceElectrons(valElectrons), d_sElectrons(sElectrons),
                        d_pElectrons(pElectrons), d_dElectrons(dElectrons),d_fElectrons(fElectrons),
                        d_color(color) {};

                /* Element symbol as string */
                std::string symbol() const { return d_symbol; }

                /* Atomic number */
                unsigned int Z() const { return d_Z; }

                /* Mass in atomic mass units (u) */
                double mass() const { return d_mass; }

                /* Van-der-Waals Radius in atomic units */
                double vdWRadius() const {
                    if (d_vdWRadius > 0.0)
                        return d_vdWRadius;
                    else throw DataNotAvailable();
                }

                /* Number of valence electrons */
                int valElectrons() const {
                    if (d_valenceElectrons > -1)
                        return d_valenceElectrons;
                    else throw DataNotAvailable();
                }

                int innerShellElectrons() const {
                    if (d_valenceElectrons > -1)
                        return Z()-d_valenceElectrons;
                    else throw DataNotAvailable();
                }

                /* Number of s-valence electrons */
                int sElectrons() const {
                    if (d_sElectrons > -1)
                        return d_sElectrons;
                    else throw DataNotAvailable();
                }

                /* Number of p-valence electrons */
                int pElectrons() const {
                    if (d_pElectrons > -1)
                        return d_pElectrons;
                    else throw DataNotAvailable();
                }

                /* Number of d-valence electrons */
                int dElectrons() const {
                    if (d_dElectrons > -1)
                        return d_dElectrons;
                    else throw DataNotAvailable();
                }

                /* Number of d-valence electrons */
                int fElectrons() const {
                    if (d_fElectrons > -1)
                        return d_fElectrons;
                    else throw DataNotAvailable();
                }

                /* Element ElememtColor for plotting */
                ElementColor color() const {
                    if (d_color.R != 0 && d_color.B != 0 && d_color.B != 0)
                        return d_color;
                    else throw DataNotAvailable();
                }

            private:
                std::string d_symbol;

                unsigned int d_Z;
                double d_mass;
                double d_vdWRadius;

                int d_valenceElectrons;
                int d_sElectrons;
                int d_pElectrons;
                int d_dElectrons;
                int d_fElectrons;

                ElementColor d_color;
            };

        private:
            /* Types internally used. */
            typedef std::map<ElementType, ElementData> ContainerType;
            typedef ContainerType::value_type ValueType;

        public:
            /* Exception for range check */
            class ElementSymbolNotFound {};

            /* Access element information based on type. Fastest lookup.
               Range checked. Throws std::out_of_range exception if no entry is available. */
            const ElementData& operator[](const ElementType& type) const;

            /* Access element information based on element symbol. Slow lookup.
               Returns a std::pair<ElementType, ElementData> or throws ElementSymbolNotFound exception. */
            const ValueType& operator[](const std::string& symbol) const;

            /* Access singleton instance. Creates one if necessary. */
            static const pElementInfo& instance();

        private:
            /* Private constructor to prevent instantiation by client. */
            pElementInfo();

            /* Keep singleton instance */
            static pElementInfo* d_instance;

            /* Store internal map */
            ContainerType d_container;

            /* create hard coded element data */
            void init_data();

        };

    } // namespace internal

} // namespace Elements

#endif // INPSIGHTS_PELEMENTINFO_H
