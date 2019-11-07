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

#include <IntegerSortedTreeWidgetItem.h>

#include "IntegerSortedTreeWidgetItem.h"

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem()
        : QTreeWidgetItem() {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidget *tree)
        : QTreeWidgetItem(tree) {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidget *parent, const QStringList &strings)
        : QTreeWidgetItem(parent, strings) {}

IntegerSortedTreeWidgetItem::IntegerSortedTreeWidgetItem(QTreeWidgetItem *parent, const QStringList &strings, int type)
        : QTreeWidgetItem(parent, strings, type) {}

bool IntegerSortedTreeWidgetItem::operator<(const QTreeWidgetItem &other) const {
    auto sortCol = treeWidget()->sortColumn();
    auto myNumber = text(sortCol).toDouble(); //TODO Refactor
    auto  otherNumber = other.text(sortCol).toDouble();
    return myNumber < otherNumber;
}