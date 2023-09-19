# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/coevolution/visual/CoevolutionViewer.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_CoevolutionViewer(object):
    def setupUi(self, CoevolutionViewer):
        CoevolutionViewer.setObjectName("CoevolutionViewer")
        CoevolutionViewer.resize(538, 444)
        self.verticalLayout = QtWidgets.QVBoxLayout(CoevolutionViewer)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label_3 = QtWidgets.QLabel(CoevolutionViewer)
        self.label_3.setObjectName("label_3")
        self.verticalLayout.addWidget(self.label_3)
        self.horizontalLayout = QtWidgets.QHBoxLayout()
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.label = QtWidgets.QLabel(CoevolutionViewer)
        self.label.setObjectName("label")
        self.horizontalLayout.addWidget(self.label)
        self.structureCombo = QtWidgets.QComboBox(CoevolutionViewer)
        self.structureCombo.setObjectName("structureCombo")
        self.horizontalLayout.addWidget(self.structureCombo)
        self.refreshButton = QtWidgets.QPushButton(CoevolutionViewer)
        self.refreshButton.setObjectName("refreshButton")
        self.horizontalLayout.addWidget(self.refreshButton)
        self.selectionErrorLabel = QtWidgets.QLabel(CoevolutionViewer)
        self.selectionErrorLabel.setStyleSheet("color: red;")
        self.selectionErrorLabel.setText("")
        self.selectionErrorLabel.setObjectName("selectionErrorLabel")
        self.horizontalLayout.addWidget(self.selectionErrorLabel)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout.addItem(spacerItem)
        self.verticalLayout.addLayout(self.horizontalLayout)
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout()
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.label_2 = QtWidgets.QLabel(CoevolutionViewer)
        self.label_2.setObjectName("label_2")
        self.horizontalLayout_2.addWidget(self.label_2)
        self.sequenceCombo = QtWidgets.QComboBox(CoevolutionViewer)
        self.sequenceCombo.setObjectName("sequenceCombo")
        self.horizontalLayout_2.addWidget(self.sequenceCombo)
        spacerItem1 = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.horizontalLayout_2.addItem(spacerItem1)
        self.verticalLayout.addLayout(self.horizontalLayout_2)
        self.tabWidget = QtWidgets.QTabWidget(CoevolutionViewer)
        self.tabWidget.setObjectName("tabWidget")
        self.verticalLayout.addWidget(self.tabWidget)

        self.retranslateUi(CoevolutionViewer)
        self.tabWidget.setCurrentIndex(-1)
        QtCore.QMetaObject.connectSlotsByName(CoevolutionViewer)

    def retranslateUi(self, CoevolutionViewer):
        _translate = QtCore.QCoreApplication.translate
        CoevolutionViewer.setWindowTitle(_translate("CoevolutionViewer", "Form"))
        self.label_3.setText(_translate("CoevolutionViewer", "Select visualization structure"))
        self.label.setText(_translate("CoevolutionViewer", "structure"))
        self.refreshButton.setText(_translate("CoevolutionViewer", "Refresh"))
        self.label_2.setText(_translate("CoevolutionViewer", "structure sequence"))
