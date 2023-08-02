# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/chimeras/visual/ChimerasGenerator.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_ChimerasGenerator(object):
    def setupUi(self, ChimerasGenerator):
        ChimerasGenerator.setObjectName("ChimerasGenerator")
        ChimerasGenerator.resize(604, 518)
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(ChimerasGenerator)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout()
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label_2 = QtWidgets.QLabel(ChimerasGenerator)
        font = QtGui.QFont()
        font.setBold(True)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setObjectName("verticalLayout")
        self.fragmentsList = QtWidgets.QListWidget(ChimerasGenerator)
        self.fragmentsList.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.fragmentsList.setObjectName("fragmentsList")
        self.verticalLayout.addWidget(self.fragmentsList)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem)
        self.generateButton = QtWidgets.QPushButton(ChimerasGenerator)
        self.generateButton.setObjectName("generateButton")
        self.horizontalLayout_2.addWidget(self.generateButton)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.fragmentsStack = QtWidgets.QStackedWidget(ChimerasGenerator)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.fragmentsStack.sizePolicy().hasHeightForWidth())
        self.fragmentsStack.setSizePolicy(sizePolicy)
        self.fragmentsStack.setObjectName("fragmentsStack")
        self.horizontalLayout.addWidget(self.fragmentsStack)
        self.horizontalLayout.setStretch(0, 1)
        self.horizontalLayout.setStretch(1, 4)
        self.verticalLayout_2.addLayout(self.horizontalLayout)
        self.verticalLayout_4.addLayout(self.verticalLayout_2)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label = QtWidgets.QLabel(ChimerasGenerator)
        font = QtGui.QFont()
        font.setBold(True)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout_3.addWidget(self.label)
        self.resultsTable = QtWidgets.QTableWidget(ChimerasGenerator)
        self.resultsTable.setObjectName("resultsTable")
        self.resultsTable.setColumnCount(1)
        self.resultsTable.setRowCount(0)
        item = QtWidgets.QTableWidgetItem()
        self.resultsTable.setHorizontalHeaderItem(0, item)
        self.verticalLayout_3.addWidget(self.resultsTable)
        self.verticalLayout_4.addLayout(self.verticalLayout_3)
        self.verticalLayout_4.setStretch(0, 1)
        self.verticalLayout_4.setStretch(1, 2)

        self.retranslateUi(ChimerasGenerator)
        QtCore.QMetaObject.connectSlotsByName(ChimerasGenerator)

    def retranslateUi(self, ChimerasGenerator):
        _translate = QtCore.QCoreApplication.translate
        ChimerasGenerator.setWindowTitle(_translate("ChimerasGenerator", "Form"))
        self.label_2.setText(_translate("ChimerasGenerator", "Fragments"))
        self.generateButton.setText(_translate("ChimerasGenerator", "Generate"))
        self.label.setText(_translate("ChimerasGenerator", "Results"))
        item = self.resultsTable.horizontalHeaderItem(0)
        item.setText(_translate("ChimerasGenerator", "Sequence"))