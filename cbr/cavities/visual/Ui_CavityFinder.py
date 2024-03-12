# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/cavities/visual/CavityFinder.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_CavityFinder(object):
    def setupUi(self, CavityFinder):
        CavityFinder.setObjectName("CavityFinder")
        CavityFinder.resize(891, 563)
        self.gridLayout = QtWidgets.QGridLayout(CavityFinder)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(CavityFinder)
        self.label.setStyleSheet("font-size: 26px;\n"
"fond-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 4)
        self.label_3 = QtWidgets.QLabel(CavityFinder)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 1, 4, 1, 1)
        self.maxSizeSelector = QtWidgets.QSpinBox(CavityFinder)
        self.maxSizeSelector.setMinimum(1)
        self.maxSizeSelector.setMaximum(999999)
        self.maxSizeSelector.setProperty("value", 125)
        self.maxSizeSelector.setObjectName("maxSizeSelector")
        self.gridLayout.addWidget(self.maxSizeSelector, 1, 7, 1, 1)
        self.busyProgress = QtWidgets.QProgressBar(CavityFinder)
        self.busyProgress.setMaximum(0)
        self.busyProgress.setProperty("value", -1)
        self.busyProgress.setObjectName("busyProgress")
        self.gridLayout.addWidget(self.busyProgress, 3, 0, 1, 15)
        self.label_4 = QtWidgets.QLabel(CavityFinder)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 1, 6, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(228, 21, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 1, 11, 1, 4)
        self.minSizeSelector = QtWidgets.QSpinBox(CavityFinder)
        self.minSizeSelector.setMinimum(0)
        self.minSizeSelector.setMaximum(999)
        self.minSizeSelector.setObjectName("minSizeSelector")
        self.gridLayout.addWidget(self.minSizeSelector, 1, 5, 1, 1)
        self.findButton = QtWidgets.QPushButton(CavityFinder)
        self.findButton.setObjectName("findButton")
        self.gridLayout.addWidget(self.findButton, 1, 9, 1, 1)
        self.refreshButton = QtWidgets.QPushButton(CavityFinder)
        self.refreshButton.setObjectName("refreshButton")
        self.gridLayout.addWidget(self.refreshButton, 1, 3, 1, 1)
        self.label_2 = QtWidgets.QLabel(CavityFinder)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 0, 1, 1)
        self.structureCombo = QtWidgets.QComboBox(CavityFinder)
        self.structureCombo.setMinimumSize(QtCore.QSize(150, 0))
        self.structureCombo.setObjectName("structureCombo")
        self.gridLayout.addWidget(self.structureCombo, 1, 1, 1, 2)
        self.cavitiesTable = QtWidgets.QTableView(CavityFinder)
        self.cavitiesTable.setObjectName("cavitiesTable")
        self.gridLayout.addWidget(self.cavitiesTable, 2, 0, 1, 15)
        self.advancedSettingsButton = QtWidgets.QPushButton(CavityFinder)
        self.advancedSettingsButton.setObjectName("advancedSettingsButton")
        self.gridLayout.addWidget(self.advancedSettingsButton, 1, 10, 1, 1)

        self.retranslateUi(CavityFinder)
        QtCore.QMetaObject.connectSlotsByName(CavityFinder)

    def retranslateUi(self, CavityFinder):
        _translate = QtCore.QCoreApplication.translate
        CavityFinder.setWindowTitle(_translate("CavityFinder", "Form"))
        self.label.setText(_translate("CavityFinder", "Cavity Finder"))
        self.label_3.setText(_translate("CavityFinder", "min size (Å)"))
        self.label_4.setText(_translate("CavityFinder", "max size (Å)"))
        self.findButton.setText(_translate("CavityFinder", "Find"))
        self.refreshButton.setText(_translate("CavityFinder", "Refresh"))
        self.label_2.setText(_translate("CavityFinder", "Structure"))
        self.advancedSettingsButton.setText(_translate("CavityFinder", "Advanced Options"))
