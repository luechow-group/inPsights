/* Copyright (C) 2018-2019 Michael Heuer.
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

#include <InPsightsWidget.h>
#include <Version.h>
#include <QApplication>
#include <memory>
#include <spdlog/spdlog.h>

int main(int argc, char *argv[]) {

    Q_INIT_RESOURCE(myresources);
    QApplication app(argc, argv);
    setlocale(LC_NUMERIC,"C");

    app.setWindowIcon(QIcon(":inPsightsIcon.png"));
    app.setApplicationName("inPsights");

    spdlog::info("Welcome to inPsights (Version: {})!", inPsights::version);

    std::unique_ptr<InPsightsWidget> widget;

    if (argc == 2)
        widget = std::make_unique<InPsightsWidget>(nullptr, argv[1]);
    else
        widget = std::make_unique<InPsightsWidget>();

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
