#include "delete.h"

MsaSelectStructure::MsaSelectStructure(QObject *parent)
    : QAbstractTableModel(parent)
{
}

QVariant MsaSelectStructure::headerData(int section, Qt::Orientation orientation, int role) const
{
    // FIXME: Implement me!
}

int MsaSelectStructure::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    // FIXME: Implement me!
}

int MsaSelectStructure::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    // FIXME: Implement me!
}

QVariant MsaSelectStructure::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    // FIXME: Implement me!
    return QVariant();
}
