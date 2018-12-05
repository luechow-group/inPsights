//
// Created by heuer on 03.12.18.
//

#ifndef INPSIGHTS_ICONNECTION_H
#define INPSIGHTS_ICONNECTION_H

#include <Qt3DCore/QEntity>

class IConnection : public Qt3DCore::QEntity {
public:
    explicit IConnection(Qt3DCore::QEntity *parent = nullptr)
            : Qt3DCore::QEntity(parent) {}
};


#endif //INPSIGHTS_ICONNECTION_H
