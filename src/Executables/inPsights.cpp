//
// Created by heuer on 15.11.18.
//

#include <InPsightsWidget.h>
#include <QApplication>

int main(int argc, char *argv[]) {

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    new InPsightsWidget();

    return QApplication::exec();

    /* TODO
     * METHOD
     * - spatial permutations (by value range or struct sim)
     * BASELIB
     * - debug linkedParticles
     * GUI
     * - state machine for buttons
     * - show correlation lines
     * - tree view
     * - energy view
     * - mouse events
     *  - energy view triggers atom hightlight
     * - dotted and dashed lines
     * - founds + psi
     * - show energies
     * - show eigenvectors (.wf needed)
     * - load screen
     * - screenshot button
     * - camera reset and settings in MoleculeWidget
     */
}
