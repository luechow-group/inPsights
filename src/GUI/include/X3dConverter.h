//
// Created by Michael Heuer on 07.10.19.
//


#ifndef INPSIGHTS_X3DCONVERTER_H
#define INPSIGHTS_X3DCONVERTER_H

#include <string>
#include <iostream>
#include <fstream>
#include "Sphere.h"
#include "Cylinder.h"

class X3dConverter{
public:
    X3dConverter(const std::string &filename,
                 const std::string& title = "Title",
                 const std::string& comment = "Comment");

    void startScene(const std::string& title, const std::string& comment);

    void addSphere(const Sphere & sphere, unsigned sortKey = 1);
    void addCylinder(const Cylinder & cylinder, unsigned sortKey = 1);
    void closeScene();

private:
    std::ofstream file;
};

#endif //INPSIGHTS_X3DCONVERTER_H
