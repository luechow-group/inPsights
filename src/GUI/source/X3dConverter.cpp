//
// Created by heuer on 08.10.19.
//

#include <X3dConverter.h>

X3dConverter::X3dConverter(const std::string &filename,
             const std::string& title,
             const std::string& comment) {
    file.open(filename);
    startScene(title, comment);
}

void X3dConverter::startScene(const std::string& title, const std::string& comment){
    file << "<html> \n"
            "<head> \n"
            "<title>" + title + "</title>\n"
                                "<script type='text/javascript' src='http://www.x3dom.org/download/x3dom.js'> </script>\n"
                                "<link rel='stylesheet' type='text/css' href='http://www.x3dom.org/download/x3dom.css'></link>\n"
                                "</head>\n"
                                "<body>\n"
                                "<h1>" + title + "</h1>\n"
                                                 "<p>" + comment + "</p>\n"
                                                                   "<center> \n"
                                                                   "<x3d width='600px' height='400px'> \n"
                                                                   "<scene>\n"
                                                                   "<DirectionalLight direction='0,0,-1' global='true'></DirectionalLight>\n";
};

void X3dConverter::addSphere(const Sphere & sphere, unsigned sortKey){
    sphere.addToXml(file,sortKey);
}

void X3dConverter::addCylinder(const Cylinder & cylinder, unsigned sortKey){
    cylinder.addToXml(file, sortKey);
}

void X3dConverter::closeScene(){
    file << "</scene>\n"
            "</x3d>\n"
            "<center>\n"
            "</body>\n"
            "</html>";
    file.close();
};