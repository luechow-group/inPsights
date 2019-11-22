/* Copyright (C) 2019 Michael Heuer.
 *
 * This file is part of inPsights.
 * inPsights is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * inPsights is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with inPsights. If not, see <https://www.gnu.org/licenses/>.
 */

#include <X3domConverter.h>


X3domConverter::X3domConverter(const std::string &filename,
                               const std::string& title,
                               const std::string& comment) {
    file.open(filename);
    startScene(title, comment);
}
X3domConverter::X3domConverter(const std::string &filename)
    : X3domConverter(filename, filename, ""){}


void X3domConverter::startScene(const std::string& title, const std::string& comment){
    file << "<html> \n"
            "<head> \n"
            "<title>" + title + "</title>\n"
            "<script type='text/javascript' src='http://www.x3dom.org/download/x3dom.js'> </script>\n"
            "<link rel='stylesheet' type='text/css' href='http://www.x3dom.org/download/x3dom.css'></link>\n"
            "</head>\n"
            "<body>\n"
            /*"<h1>" + title + "</h1>\n"
            "<p>" + comment + "</p>\n"*/
            "<center> \n"
            "<x3d width='100%' height='100%'> \n"
            "<scene>\n"
            "<DirectionalLight direction='0,0,-1' global='true'></DirectionalLight>\n";
};

void X3domConverter::addSphere(const Sphere & sphere, unsigned sortKey){
    sphere.addToXml(file,sortKey);
}

void X3domConverter::addCylinder(const Cylinder & cylinder, unsigned sortKey){
    cylinder.addToXml(file, sortKey);
}

void X3domConverter::closeScene(){
    file << "</scene>\n"
            "</x3d>\n"
            "</center>\n"
            "</body>\n"
            "</html>";
    file.close();
};