// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef INPSIGHTS_ICONNECTION_H
#define INPSIGHTS_ICONNECTION_H

#include <Qt3DCore/QEntity>

class IConnection : public Qt3DCore::QEntity {
public:
    explicit IConnection(Qt3DCore::QEntity *parent = nullptr)
            : Qt3DCore::QEntity(parent) {}
};


#endif //INPSIGHTS_ICONNECTION_H
