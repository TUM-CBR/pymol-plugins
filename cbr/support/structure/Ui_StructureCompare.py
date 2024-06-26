# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/support/structure/StructureCompare.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_StructureCompare(object):
    def setupUi(self, StructureCompare):
        StructureCompare.setObjectName("StructureCompare")
        StructureCompare.resize(579, 464)
        self.gridLayout = QtWidgets.QGridLayout(StructureCompare)
        self.gridLayout.setObjectName("gridLayout")
        self.structuresTable = QtWidgets.QTableView(StructureCompare)
        self.structuresTable.setObjectName("structuresTable")
        self.gridLayout.addWidget(self.structuresTable, 1, 0, 1, 4)
        self.structureColoringCombo = QtWidgets.QComboBox(StructureCompare)
        self.structureColoringCombo.setObjectName("structureColoringCombo")
        self.gridLayout.addWidget(self.structureColoringCombo, 2, 2, 1, 1)
        self.selectStructuresButton = QtWidgets.QPushButton(StructureCompare)
        self.selectStructuresButton.setObjectName("selectStructuresButton")
        self.gridLayout.addWidget(self.selectStructuresButton, 2, 3, 1, 1)
        self.label_2 = QtWidgets.QLabel(StructureCompare)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 1, 1, 1)
        self.label = QtWidgets.QLabel(StructureCompare)
        self.label.setStyleSheet("font-size: 24pt;\n"
"font-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 4)
        self.gridLayout.setColumnStretch(0, 1)

        self.retranslateUi(StructureCompare)
        QtCore.QMetaObject.connectSlotsByName(StructureCompare)

    def retranslateUi(self, StructureCompare):
        _translate = QtCore.QCoreApplication.translate
        StructureCompare.setWindowTitle(_translate("StructureCompare", "Form"))
        self.selectStructuresButton.setText(_translate("StructureCompare", "Select Structures"))
        self.label_2.setText(_translate("StructureCompare", "Structure Coloring"))
        self.label.setText(_translate("StructureCompare", "Structure Compare"))
