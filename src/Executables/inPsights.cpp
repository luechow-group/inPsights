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
     * - investigate maxima sith two-electrons at hydrogen core (gradient might be two small and 
     * electrons slightly too far away from nucleus)
     * - spatial permutations (by value range or struct sim)
     * - investigate N2H4 hydrogen: two electrons close at hydrogen, 1 at core, other one nearby
     *
     * BASELIB
     * - Molecule class (atomsVector, electronsVectorCollection)
     * -
     * GUI
     *  - dotted and dashed lines
     * - select clusters (list with checkboxes, max selected structures 100?)
     * - founds + psi
     * - show energies
     * - show eigenvectors (.wf needed)
     */
}
