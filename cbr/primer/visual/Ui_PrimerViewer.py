# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/primer/visual/PrimerViewer.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_PrimerViewer(object):
    def setupUi(self, PrimerViewer):
        PrimerViewer.setObjectName("PrimerViewer")
        PrimerViewer.resize(650, 622)
        self.verticalLayout = QtWidgets.QVBoxLayout(PrimerViewer)
        self.verticalLayout.setObjectName("verticalLayout")
        self.sequenceTable = QtWidgets.QTableWidget(PrimerViewer)
        self.sequenceTable.setObjectName("sequenceTable")
        self.sequenceTable.setColumnCount(0)
        self.sequenceTable.setRowCount(0)
        self.verticalLayout.addWidget(self.sequenceTable)
        self.primersTable = QtWidgets.QTableWidget(PrimerViewer)
        self.primersTable.setObjectName("primersTable")
        self.primersTable.setColumnCount(0)
        self.primersTable.setRowCount(0)
        self.verticalLayout.addWidget(self.primersTable)
        self.verticalLayout.setStretch(0, 1)

        self.retranslateUi(PrimerViewer)
        QtCore.QMetaObject.connectSlotsByName(PrimerViewer)

    def retranslateUi(self, PrimerViewer):
        _translate = QtCore.QCoreApplication.translate
        PrimerViewer.setWindowTitle(_translate("PrimerViewer", "Form"))
