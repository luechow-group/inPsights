//
// Created by Michael Heuer on 17.11.18.
//

#ifndef INPSIGHTS_INPSIGHTSSPLASHSCREEN_H
#define INPSIGHTS_INPSIGHTSSPLASHSCREEN_H

#include <QSplashScreen>

namespace InPsightsSplashScreen {
    QSplashScreen* getInPsightsSplashScreen(){
        QSplashScreen * splash = new QSplashScreen;

        QImage file(":inPsights.png");
        QPixmap qPixmap;
        qPixmap.convertFromImage(file);
        splash->setPixmap(qPixmap);
        splash->show();
        splash->showMessage("Version 1.0.0\nby M. Heuer", Qt::AlignRight, Qt::gray);

        return splash;
    }
}

#endif //INPSIGHTS_INPSIGHTSSPLASHSCREEN_H
