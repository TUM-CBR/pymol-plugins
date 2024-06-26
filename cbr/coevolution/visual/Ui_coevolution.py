# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/coevolution/visual/coevolution.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_Coevolution(object):
    def setupUi(self, Coevolution):
        Coevolution.setObjectName("Coevolution")
        Coevolution.resize(987, 612)
        self.gridLayout = QtWidgets.QGridLayout(Coevolution)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(Coevolution)
        self.label.setStyleSheet("font-size: 26px;\n"
"font-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 2)
        self.selectAlignmentButton = QtWidgets.QPushButton(Coevolution)
        self.selectAlignmentButton.setObjectName("selectAlignmentButton")
        self.gridLayout.addWidget(self.selectAlignmentButton, 1, 0, 1, 1)
        self.selectedAlignmentLabel = QtWidgets.QLabel(Coevolution)
        self.selectedAlignmentLabel.setMinimumSize(QtCore.QSize(100, 0))
        self.selectedAlignmentLabel.setText("")
        self.selectedAlignmentLabel.setObjectName("selectedAlignmentLabel")
        self.gridLayout.addWidget(self.selectedAlignmentLabel, 1, 1, 1, 1)
        self.selectStructureButton = QtWidgets.QPushButton(Coevolution)
        self.selectStructureButton.setObjectName("selectStructureButton")
        self.gridLayout.addWidget(self.selectStructureButton, 1, 2, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout.addItem(spacerItem, 1, 3, 1, 1)
        self.parametersTable = QtWidgets.QTableView(Coevolution)
        self.parametersTable.setMaximumSize(QtCore.QSize(16777215, 60))
        self.parametersTable.setObjectName("parametersTable")
        self.gridLayout.addWidget(self.parametersTable, 1, 4, 2, 1)
        self.label_2 = QtWidgets.QLabel(Coevolution)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 0, 1, 1)
        self.resultsNumberBox = QtWidgets.QSpinBox(Coevolution)
        self.resultsNumberBox.setMinimumSize(QtCore.QSize(0, 0))
        self.resultsNumberBox.setMinimum(1)
        self.resultsNumberBox.setMaximum(500)
        self.resultsNumberBox.setProperty("value", 50)
        self.resultsNumberBox.setObjectName("resultsNumberBox")
        self.gridLayout.addWidget(self.resultsNumberBox, 2, 1, 1, 1)
        self.onlyStructureCombo = QtWidgets.QCheckBox(Coevolution)
        self.onlyStructureCombo.setObjectName("onlyStructureCombo")
        self.gridLayout.addWidget(self.onlyStructureCombo, 2, 2, 1, 1)
        self.splitter = QtWidgets.QSplitter(Coevolution)
        self.splitter.setOrientation(QtCore.Qt.Orientation.Horizontal)
        self.splitter.setObjectName("splitter")
        self.alignmentTable = QtWidgets.QTableView(self.splitter)
        self.alignmentTable.setObjectName("alignmentTable")
        self.detailsTable = QtWidgets.QTableView(self.splitter)
        self.detailsTable.setObjectName("detailsTable")
        self.gridLayout.addWidget(self.splitter, 3, 0, 1, 5)
        self.busyProgress = QtWidgets.QProgressBar(Coevolution)
        self.busyProgress.setMaximum(0)
        self.busyProgress.setProperty("value", -1)
        self.busyProgress.setObjectName("busyProgress")
        self.gridLayout.addWidget(self.busyProgress, 4, 0, 1, 5)
        self.gridLayout.setColumnStretch(3, 1)
        self.gridLayout.setColumnStretch(4, 2)
        self.gridLayout.setRowStretch(3, 1)

        self.retranslateUi(Coevolution)
        QtCore.QMetaObject.connectSlotsByName(Coevolution)

    def retranslateUi(self, Coevolution):
        _translate = QtCore.QCoreApplication.translate
        Coevolution.setWindowTitle(_translate("Coevolution", "Form"))
        self.label.setText(_translate("Coevolution", "Coevolution"))
        self.selectAlignmentButton.setText(_translate("Coevolution", "Select Alignment"))
        self.selectStructureButton.setText(_translate("Coevolution", "Select Structures"))
        self.label_2.setText(_translate("Coevolution", "Results number"))
        self.onlyStructureCombo.setText(_translate("Coevolution", "Scope Structure"))
