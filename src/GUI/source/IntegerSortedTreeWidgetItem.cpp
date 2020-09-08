// Copyright (C) 2018-2019 Michael Heuer.
// SPDX-License-Identifier: GPL-3.0-or-later

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