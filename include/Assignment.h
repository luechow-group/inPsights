//
// Created by Morian Sonnet on 19.05.2017.
//

#ifndef LOCALSPINMULTIPLICITY_ASSIGNMENT_H
#define LOCALSPINMULTIPLICITY_ASSIGNMENT_H
#include <vector>
#include <cstdlib>

/*
 * Data type for holding an Electron-Core-Assignment
 * The top level vector contains of the Electron-Core-Assignment for a single Core.
 * The Electron-Core-Assignment for a single Core is a pair of the specific Core and a Vector of Electrons assigned to it.
 */
typedef std::vector<std::pair<int,std::vector<int> > >  Assignment;

#endif //LOCALSPINMULTIPLICITY_ASSIGNMENT_H
