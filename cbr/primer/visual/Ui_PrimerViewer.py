# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/primer/visual/PrimerViewer.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_PrimerViewer(object):
    def setupUi(self, PrimerViewer):
        PrimerViewer.setObjectName("PrimerViewer")
        PrimerViewer.resize(877, 733)
        self.verticalLayout = QtWidgets.QVBoxLayout(PrimerViewer)
        self.verticalLayout.setObjectName("verticalLayout")
        self.resultsWidget = QtWidgets.QWidget(PrimerViewer)
        self.resultsWidget.setObjectName("resultsWidget")
        self.verticalLayout.addWidget(self.resultsWidget)
        self.primersTable = QtWidgets.QTableView(PrimerViewer)
        self.primersTable.setObjectName("primersTable")
        self.verticalLayout.addWidget(self.primersTable)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.tmWeightEdit = QtWidgets.QLineEdit(PrimerViewer)
        self.tmWeightEdit.setObjectName("tmWeightEdit")
        self.gridLayout.addWidget(self.tmWeightEdit, 5, 1, 1, 1)
        self.tmPrimersWeightEdit = QtWidgets.QLineEdit(PrimerViewer)
        self.tmPrimersWeightEdit.setObjectName("tmPrimersWeightEdit")
        self.gridLayout.addWidget(self.tmPrimersWeightEdit, 5, 5, 1, 1)
        self.label_6 = QtWidgets.QLabel(PrimerViewer)
        self.label_6.setObjectName("label_6")
        self.gridLayout.addWidget(self.label_6, 5, 4, 1, 1)
        self.label_5 = QtWidgets.QLabel(PrimerViewer)
        self.label_5.setObjectName("label_5")
        self.gridLayout.addWidget(self.label_5, 4, 4, 1, 1)
        self.label = QtWidgets.QLabel(PrimerViewer)
        self.label.setStyleSheet("font-weight: bold")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 2, 0, 1, 1)
        self.tmDeltaWeightEdit = QtWidgets.QLineEdit(PrimerViewer)
        self.tmDeltaWeightEdit.setObjectName("tmDeltaWeightEdit")
        self.gridLayout.addWidget(self.tmDeltaWeightEdit, 4, 8, 1, 1)
        self.label_2 = QtWidgets.QLabel(PrimerViewer)
        self.label_2.setStyleSheet("font-weight: bold;")
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 2, 4, 1, 1)
        self.tmEdit = QtWidgets.QLineEdit(PrimerViewer)
        self.tmEdit.setObjectName("tmEdit")
        self.gridLayout.addWidget(self.tmEdit, 4, 1, 1, 2)
        self.label_3 = QtWidgets.QLabel(PrimerViewer)
        self.label_3.setObjectName("label_3")
        self.gridLayout.addWidget(self.label_3, 5, 0, 1, 1)
        self.selectButton = QtWidgets.QPushButton(PrimerViewer)
        self.selectButton.setObjectName("selectButton")
        self.gridLayout.addWidget(self.selectButton, 5, 9, 1, 1)
        self.label_7 = QtWidgets.QLabel(PrimerViewer)
        self.label_7.setStyleSheet("font-weight: bold;")
        self.label_7.setObjectName("label_7")
        self.gridLayout.addWidget(self.label_7, 2, 7, 1, 1)
        self.tmPrimersEdit = QtWidgets.QLineEdit(PrimerViewer)
        self.tmPrimersEdit.setObjectName("tmPrimersEdit")
        self.gridLayout.addWidget(self.tmPrimersEdit, 4, 5, 1, 1)
        self.exportTableButton = QtWidgets.QPushButton(PrimerViewer)
        self.exportTableButton.setObjectName("exportTableButton")
        self.gridLayout.addWidget(self.exportTableButton, 4, 9, 1, 1)
        self.label_8 = QtWidgets.QLabel(PrimerViewer)
        self.label_8.setObjectName("label_8")
        self.gridLayout.addWidget(self.label_8, 4, 7, 1, 1)
        self.label_4 = QtWidgets.QLabel(PrimerViewer)
        self.label_4.setObjectName("label_4")
        self.gridLayout.addWidget(self.label_4, 4, 0, 1, 1)
        self.saveAllButton = QtWidgets.QPushButton(PrimerViewer)
        self.saveAllButton.setObjectName("saveAllButton")
        self.gridLayout.addWidget(self.saveAllButton, 2, 9, 1, 1)
        self.filterProgress = QtWidgets.QProgressBar(PrimerViewer)
        self.filterProgress.setMinimumSize(QtCore.QSize(0, 25))
        self.filterProgress.setStyleSheet("")
        self.filterProgress.setProperty("value", 24)
        self.filterProgress.setObjectName("filterProgress")
        self.gridLayout.addWidget(self.filterProgress, 6, 1, 1, 9)
        self.widget = QtWidgets.QWidget(PrimerViewer)
        self.widget.setMinimumSize(QtCore.QSize(0, 25))
        self.widget.setObjectName("widget")
        self.gridLayout.addWidget(self.widget, 6, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout)
        self.verticalLayout.setStretch(0, 5)

        self.retranslateUi(PrimerViewer)
        QtCore.QMetaObject.connectSlotsByName(PrimerViewer)

    def retranslateUi(self, PrimerViewer):
        _translate = QtCore.QCoreApplication.translate
        PrimerViewer.setWindowTitle(_translate("PrimerViewer", "Form"))
        self.label_6.setText(_translate("PrimerViewer", "Weight"))
        self.label_5.setText(_translate("PrimerViewer", "Value"))
        self.label.setText(_translate("PrimerViewer", "Tm (all)"))
        self.label_2.setText(_translate("PrimerViewer", "Tm (primers)"))
        self.label_3.setText(_translate("PrimerViewer", "Weight"))
        self.selectButton.setText(_translate("PrimerViewer", "Select"))
        self.label_7.setText(_translate("PrimerViewer", "Tm (delta)"))
        self.exportTableButton.setText(_translate("PrimerViewer", "Export Current Results"))
        self.label_8.setText(_translate("PrimerViewer", "Weight"))
        self.label_4.setText(_translate("PrimerViewer", "Value"))
        self.saveAllButton.setText(_translate("PrimerViewer", "Save All Primers"))
