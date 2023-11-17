# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/support/msa/visual/MsaStructureSelector.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MsaStructureSelector(object):
    def setupUi(self, MsaStructureSelector):
        MsaStructureSelector.setObjectName("MsaStructureSelector")
        MsaStructureSelector.resize(518, 138)
        self.gridLayout = QtWidgets.QGridLayout(MsaStructureSelector)
        self.gridLayout.setObjectName("gridLayout")
        self.refreshButton = QtWidgets.QPushButton(MsaStructureSelector)
        self.refreshButton.setObjectName("refreshButton")
        self.gridLayout.addWidget(self.refreshButton, 1, 2, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 3, 2, 1, 1)
        self.confirmButtons = QtWidgets.QDialogButtonBox(MsaStructureSelector)
        self.confirmButtons.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.confirmButtons.setObjectName("confirmButtons")
        self.gridLayout.addWidget(self.confirmButtons, 2, 1, 1, 2)
        self.structureCombo = QtWidgets.QComboBox(MsaStructureSelector)
        self.structureCombo.setObjectName("structureCombo")
        self.gridLayout.addWidget(self.structureCombo, 1, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(325, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(MsaStructureSelector)
        self.label.setStyleSheet("font-size: 30px;\n"
"font-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 3)

        self.retranslateUi(MsaStructureSelector)
        QtCore.QMetaObject.connectSlotsByName(MsaStructureSelector)

    def retranslateUi(self, MsaStructureSelector):
        _translate = QtCore.QCoreApplication.translate
        MsaStructureSelector.setWindowTitle(_translate("MsaStructureSelector", "Dialog"))
        self.refreshButton.setText(_translate("MsaStructureSelector", "Refresh"))
        self.label.setText(_translate("MsaStructureSelector", "Select Structure of the  Sequence"))