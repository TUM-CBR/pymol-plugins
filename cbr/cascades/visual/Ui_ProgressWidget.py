# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/cascades/visual/ProgressWidget.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ProgressWidget(object):
    def setupUi(self, ProgressWidget):
        ProgressWidget.setObjectName("ProgressWidget")
        ProgressWidget.resize(541, 141)
        self.gridLayout = QtWidgets.QGridLayout(ProgressWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.progressBar = QtWidgets.QProgressBar(ProgressWidget)
        self.progressBar.setMaximum(0)
        self.progressBar.setProperty("value", -1)
        self.progressBar.setObjectName("progressBar")
        self.gridLayout.addWidget(self.progressBar, 1, 0, 1, 1)
        self.detailsButton = QtWidgets.QPushButton(ProgressWidget)
        self.detailsButton.setObjectName("detailsButton")
        self.gridLayout.addWidget(self.detailsButton, 1, 1, 1, 1)
        self.label = QtWidgets.QLabel(ProgressWidget)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.detailsTable = QtWidgets.QWidget(ProgressWidget)
        self.detailsTable.setObjectName("detailsTable")
        self.gridLayout.addWidget(self.detailsTable, 2, 0, 1, 2)

        self.retranslateUi(ProgressWidget)
        QtCore.QMetaObject.connectSlotsByName(ProgressWidget)

    def retranslateUi(self, ProgressWidget):
        _translate = QtCore.QCoreApplication.translate
        ProgressWidget.setWindowTitle(_translate("ProgressWidget", "Form"))
        self.detailsButton.setText(_translate("ProgressWidget", "Details"))
        self.label.setText(_translate("ProgressWidget", "BLASTing away, but there is a lot to BLAST. It will take some hours."))
