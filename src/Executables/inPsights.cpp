//
// Created by heuer on 15.11.18.
//

#include <InPsightsWidget.h>
#include <QApplication>
#include <Logger.h>

int main(int argc, char *argv[]) {

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    app.setWindowIcon(QIcon(":inPsightsIcon.png"));
    app.setApplicationName("inPsights");

    Logger::initialize();
    spdlog::get(Logger::name)->info("Welcome to inPsights!");


    new InPsightsWidget();

    return QApplication::exec();

    /* TODO
     *  METHOD
     *  - spatial permutations (by value range or struct sim)
     *  BASELIB
     *  - test linkedParticles
     *  GUI
     *  - state machine for buttons
     *  - mouse events
     *  - dotted and dashed lines
     *  - show eigenvectors (.wf needed)
     *  - load menu
     *  - screenshot button
     *  - axis, camera reset and settings in MoleculeWidget
     */
}
