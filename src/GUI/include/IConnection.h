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

#ifndef INPSIGHTS_ICONNECTION_H
#define INPSIGHTS_ICONNECTION_H

#include <Qt3DCore/QEntity>

class IConnection : public Qt3DCore::QEntity {
public:
    explicit IConnection(Qt3DCore::QEntity *parent = nullptr)
            : Qt3DCore::QEntity(parent) {}
};


#endif //INPSIGHTS_ICONNECTION_H
