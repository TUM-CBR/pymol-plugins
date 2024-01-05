# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/kinetics/visual/KineticsOptimizer.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_KineticsOptimizer(object):
    def setupUi(self, KineticsOptimizer):
        KineticsOptimizer.setObjectName("KineticsOptimizer")
        KineticsOptimizer.resize(1005, 610)
        self.gridLayout = QtWidgets.QGridLayout(KineticsOptimizer)
        self.gridLayout.setObjectName("gridLayout")
        self.fitProgressBar = QtWidgets.QProgressBar(KineticsOptimizer)
        self.fitProgressBar.setMaximum(0)
        self.fitProgressBar.setProperty("value", -1)
        self.fitProgressBar.setObjectName("fitProgressBar")
        self.gridLayout.addWidget(self.fitProgressBar, 7, 0, 1, 3)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 6, 2, 1, 2)
        self.parametersTable = QtWidgets.QTableView(KineticsOptimizer)
        self.parametersTable.setObjectName("parametersTable")
        self.gridLayout.addWidget(self.parametersTable, 3, 0, 2, 4)
        self.wizardButton = QtWidgets.QPushButton(KineticsOptimizer)
        self.wizardButton.setObjectName("wizardButton")
        self.gridLayout.addWidget(self.wizardButton, 6, 0, 1, 1)
        self.widget = QtWidgets.QWidget(KineticsOptimizer)
        self.widget.setMinimumSize(QtCore.QSize(0, 20))
        self.widget.setMaximumSize(QtCore.QSize(0, 16777215))
        self.widget.setObjectName("widget")
        self.gridLayout.addWidget(self.widget, 7, 3, 1, 1)
        self.plotWidget = QtWidgets.QWidget(KineticsOptimizer)
        self.plotWidget.setObjectName("plotWidget")
        self.gridLayout.addWidget(self.plotWidget, 3, 4, 6, 2)
        self.fitModelButton = QtWidgets.QPushButton(KineticsOptimizer)
        self.fitModelButton.setObjectName("fitModelButton")
        self.gridLayout.addWidget(self.fitModelButton, 6, 1, 1, 1)
        self.label = QtWidgets.QLabel(KineticsOptimizer)
        self.label.setStyleSheet("font-size: 24px;\n"
"font-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 6)
        self.label_2 = QtWidgets.QLabel(KineticsOptimizer)
        self.label_2.setStyleSheet("font-weight: bold;\n"
"font-size: 16px;")
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 1, 4, 1, 2)
        self.label_4 = QtWidgets.QLabel(KineticsOptimizer)
        self.label_4.setStyleSheet("font-weight: bold;\n"
"font-size: 16px;")
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 1, 0, 1, 4)
        self.gridLayout.setColumnStretch(5, 1)
        self.gridLayout.setRowStretch(3, 1)

        self.retranslateUi(KineticsOptimizer)
        QtCore.QMetaObject.connectSlotsByName(KineticsOptimizer)

    def retranslateUi(self, KineticsOptimizer):
        _translate = QtCore.QCoreApplication.translate
        KineticsOptimizer.setWindowTitle(_translate("KineticsOptimizer", "Form"))
        self.wizardButton.setText(_translate("KineticsOptimizer", "Parameters Wizard"))
        self.fitModelButton.setText(_translate("KineticsOptimizer", "Fit Model"))
        self.label.setText(_translate("KineticsOptimizer", "Kinetics Optimizer"))
        self.label_2.setText(_translate("KineticsOptimizer", "Fit"))
        self.label_4.setText(_translate("KineticsOptimizer", "Model"))
