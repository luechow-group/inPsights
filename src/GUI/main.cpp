//
// Created by Michael Heuer on 17.11.18.
//

#include <InPsightsSplashScreen.h>
#include <QApplication>
#include <QTimer>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC, "C");
    QWidget *window = new QWidget();
    auto splash = InPsightsSplashScreen::getInPsightsSplashScreen();

    QTimer::singleShot(4000, splash, SLOT(close()));
    QTimer::singleShot(4000, window, SLOT(show()));

    return app.exec();
}