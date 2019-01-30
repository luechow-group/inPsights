//
// Created by heuer on 15.11.18.
//

#include <InPsightsWidget.h>
#include <QApplication>
#include <memory>
#include <spdlog/spdlog.h>

int main(int argc, char *argv[]) {

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    app.setWindowIcon(QIcon(":inPsightsIcon.png"));
    app.setApplicationName("inPsights");

    spdlog::info("Welcome to inPsights!");

    auto widget = std::make_unique<InPsightsWidget>();

    return QApplication::exec();

    /* TODO
     *  METHOD
     *  - spatial permutations (by value range or struct sim)
     *  BASELIB
     *  - test linkedParticles
     *  GUI
     *  - bonds and spinCorrelations together lead to program crash (selecting correlations when bonds are displayed and vice versa)
     *  - state machine for buttons
     *  - mouse events
     *  - dotted and dashed lines
     *  - show eigenvectors (.wf needed)
     *  - load menu
     *  - screenshot button
     *  - axis, camera reset and settings in MoleculeWidget
     */
}
