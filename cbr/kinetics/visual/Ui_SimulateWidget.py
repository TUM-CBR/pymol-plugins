# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/kinetics/visual/SimulateWidget.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SimulateWidget(object):
    def setupUi(self, SimulateWidget):
        SimulateWidget.setObjectName("SimulateWidget")
        SimulateWidget.resize(982, 664)
        self.gridLayout = QtWidgets.QGridLayout(SimulateWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.fitButton = QtWidgets.QPushButton(SimulateWidget)
        self.fitButton.setObjectName("fitButton")
        self.gridLayout.addWidget(self.fitButton, 3, 2, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(151, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 3, 1, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(699, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 1, 0, 1, 1)
        self.plotWidget = QtWidgets.QWidget(SimulateWidget)
        self.plotWidget.setObjectName("plotWidget")
        self.gridLayout.addWidget(self.plotWidget, 0, 0, 1, 3)
        self.parametersWidget = QtWidgets.QTableView(SimulateWidget)
        self.parametersWidget.setObjectName("parametersWidget")
        self.gridLayout.addWidget(self.parametersWidget, 1, 1, 1, 2)
        self.gridLayout.setColumnStretch(0, 1)
        self.gridLayout.setRowStretch(0, 1)

        self.retranslateUi(SimulateWidget)
        QtCore.QMetaObject.connectSlotsByName(SimulateWidget)

    def retranslateUi(self, SimulateWidget):
        _translate = QtCore.QCoreApplication.translate
        SimulateWidget.setWindowTitle(_translate("SimulateWidget", "Form"))
        self.fitButton.setText(_translate("SimulateWidget", "Fit By Simulation"))
