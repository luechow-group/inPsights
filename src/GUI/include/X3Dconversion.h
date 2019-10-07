//
// Created by Michael Heuer on 07.10.19.
//


#ifndef INPSIGHTS_X3DCONVERSION_H
#define INPSIGHTS_X3DCONVERSION_H

#include <string>
#include <iostream>
#include <fstream>
#include "Sphere.h"
#include "Cylinder.h"

class X3Dconverter{
public:
    X3Dconverter(const std::string &filename)
    {
        file.open(filename);
        startScene();
    }

    void startScene(){
        file << "<html> \n"
                "<head> \n"
                "<title>My first X3DOM page</title>\n"
                "<script type='text/javascript' src='http://www.x3dom.org/download/x3dom.js'> </script> \n"
                "<link rel='stylesheet' type='text/css' href='http://www.x3dom.org/download/x3dom.css'></link> \n"
                "</head> \n"
                "<body> \n"
                "<h1>Hello, X3DOM!</h1> \n"
                "<p> \n"
                "This is my first html page with some 3d objects. \n"
                "</p> \n"
                "<x3d width='600px' height='400px'> \n"
                "<scene>\n";
    };

    void addSphere(const Sphere & sphere){
        file << sphere;
    }

    void addCylinder(const Cylinder & cylinder){
        file << cylinder;
    }


    void closeScene(){
        file << "</scene>\n"
                "</x3d>\n"
                "</body>\n"
                "</html>";
        file.close();
    };


private:
    std::ofstream file;

};

#endif //INPSIGHTS_X3DCONVERSION_H
