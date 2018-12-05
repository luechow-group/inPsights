//
// Created by heuer on 05.12.18.
//

#include <MaximaTreeWidgetItem.h>

#include "MaximaTreeWidgetItem.h"

MaximaTreeWidgetItem::MaximaTreeWidgetItem(QTreeWidget *tree)
        : QTreeWidgetItem(tree) {}

MaximaTreeWidgetItem::MaximaTreeWidgetItem(QTreeWidget *parent, const QStringList &strings)
        : QTreeWidgetItem(parent, strings) {}

MaximaTreeWidgetItem::MaximaTreeWidgetItem(QTreeWidgetItem *parent, const QStringList &strings, int type)
        : QTreeWidgetItem(parent, strings, type) {}

bool MaximaTreeWidgetItem::operator<(const QTreeWidgetItem &other) const {

    int sortCol = treeWidget()->sortColumn();
    int myNumber = text(sortCol).toInt();
    int otherNumber = other.text(sortCol).toInt();
    return myNumber < otherNumber;
}