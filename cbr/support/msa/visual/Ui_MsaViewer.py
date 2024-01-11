# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/support/msa/visual/MsaViewer.ui'
#
# Created by: PyQt5 UI code generator 5.15.10
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MsaViewer(object):
    def setupUi(self, MsaViewer):
        MsaViewer.setObjectName("MsaViewer")
        MsaViewer.resize(836, 538)
        self.gridLayout = QtWidgets.QGridLayout(MsaViewer)
        self.gridLayout.setObjectName("gridLayout")
        self.maskCombo = QtWidgets.QComboBox(MsaViewer)
        self.maskCombo.setMinimumSize(QtCore.QSize(150, 0))
        self.maskCombo.setMaximumSize(QtCore.QSize(150, 16777215))
        self.maskCombo.setObjectName("maskCombo")
        self.gridLayout.addWidget(self.maskCombo, 0, 1, 1, 1)
        self.label = QtWidgets.QLabel(MsaViewer)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(40, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 0, 2, 1, 1)
        self.tabWidget = QtWidgets.QTabWidget(MsaViewer)
        self.tabWidget.setObjectName("tabWidget")
        self.tableView = QtWidgets.QWidget()
        self.tableView.setObjectName("tableView")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.tableView)
        self.verticalLayout.setObjectName("verticalLayout")
        self.msaTable = QtWidgets.QTableView(self.tableView)
        self.msaTable.setObjectName("msaTable")
        self.verticalLayout.addWidget(self.msaTable)
        self.tabWidget.addTab(self.tableView, "")
        self.globalView = QtWidgets.QWidget()
        self.globalView.setObjectName("globalView")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.globalView)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.globalViewWidget = QtWidgets.QWidget(self.globalView)
        self.globalViewWidget.setObjectName("globalViewWidget")
        self.verticalLayout_2.addWidget(self.globalViewWidget)
        self.tabWidget.addTab(self.globalView, "")
        self.gridLayout.addWidget(self.tabWidget, 1, 0, 1, 3)

        self.retranslateUi(MsaViewer)
        self.tabWidget.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(MsaViewer)

    def retranslateUi(self, MsaViewer):
        _translate = QtCore.QCoreApplication.translate
        MsaViewer.setWindowTitle(_translate("MsaViewer", "Form"))
        self.label.setText(_translate("MsaViewer", "Masked Positions"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tableView), _translate("MsaViewer", "Table View"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.globalView), _translate("MsaViewer", "Global View"))
