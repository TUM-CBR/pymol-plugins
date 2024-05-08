# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/sequences/visual/SequenceSearchWidget.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SequenceSearchWidget(object):
    def setupUi(self, SequenceSearchWidget):
        SequenceSearchWidget.setObjectName("SequenceSearchWidget")
        SequenceSearchWidget.resize(882, 649)
        self.gridLayout = QtWidgets.QGridLayout(SequenceSearchWidget)
        self.gridLayout.setObjectName("gridLayout")
        self.splitter = QtWidgets.QSplitter(SequenceSearchWidget)
        self.splitter.setOrientation(QtCore.Qt.Orientation.Vertical)
        self.splitter.setObjectName("splitter")
        self.sequenceText = QtWidgets.QPlainTextEdit(self.splitter)
        self.sequenceText.setObjectName("sequenceText")
        self.resultsTable = QtWidgets.QTableView(self.splitter)
        self.resultsTable.setObjectName("resultsTable")
        self.gridLayout.addWidget(self.splitter, 1, 0, 1, 5)
        self.label = QtWidgets.QLabel(SequenceSearchWidget)
        self.label.setStyleSheet("font-size: 26px;\n"
"font-weight: bold;")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 3)
        self.searchButton = QtWidgets.QPushButton(SequenceSearchWidget)
        self.searchButton.setObjectName("searchButton")
        self.gridLayout.addWidget(self.searchButton, 2, 4, 1, 1)
        self.changeDatabase = QtWidgets.QPushButton(SequenceSearchWidget)
        self.changeDatabase.setObjectName("changeDatabase")
        self.gridLayout.addWidget(self.changeDatabase, 2, 2, 1, 1)
        self.rescanButton = QtWidgets.QPushButton(SequenceSearchWidget)
        self.rescanButton.setObjectName("rescanButton")
        self.gridLayout.addWidget(self.rescanButton, 2, 3, 1, 1)
        self.scanProgress = QtWidgets.QProgressBar(SequenceSearchWidget)
        self.scanProgress.setMaximum(0)
        self.scanProgress.setProperty("value", -1)
        self.scanProgress.setObjectName("scanProgress")
        self.gridLayout.addWidget(self.scanProgress, 2, 1, 1, 1)

        self.retranslateUi(SequenceSearchWidget)
        QtCore.QMetaObject.connectSlotsByName(SequenceSearchWidget)

    def retranslateUi(self, SequenceSearchWidget):
        _translate = QtCore.QCoreApplication.translate
        SequenceSearchWidget.setWindowTitle(_translate("SequenceSearchWidget", "Form"))
        self.label.setText(_translate("SequenceSearchWidget", "Sequence Search"))
        self.searchButton.setText(_translate("SequenceSearchWidget", "Search"))
        self.changeDatabase.setText(_translate("SequenceSearchWidget", "Change Databse"))
        self.rescanButton.setText(_translate("SequenceSearchWidget", "Rescan Files"))
