// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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

    spdlog::info("Welcome to inPsights (Version: {})!", inPsights::version());

    std::unique_ptr<InPsightsWidget> widget;

    if (argc == 2)
        widget = std::make_unique<InPsightsWidget>(nullptr, argv[1]);
    else
        widget = std::make_unique<InPsightsWidget>();

    return QApplication::exec();
}
