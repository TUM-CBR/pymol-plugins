# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'SchemaSelectorWidget.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SchemaSelectorWidget(object):
    def setupUi(self, SchemaSelectorWidget):
        SchemaSelectorWidget.setObjectName("SchemaSelectorWidget")
        SchemaSelectorWidget.resize(751, 494)
        self.verticalLayout = QtWidgets.QVBoxLayout(SchemaSelectorWidget)
        self.verticalLayout.setObjectName("verticalLayout")
        self.splitter = QtWidgets.QSplitter(SchemaSelectorWidget)
        self.splitter.setOrientation(QtCore.Qt.Horizontal)
        self.splitter.setObjectName("splitter")
        self.resultsList = QtWidgets.QListWidget(self.splitter)
        self.resultsList.setObjectName("resultsList")
        self.resultsViewer = QtWidgets.QTableWidget(self.splitter)
        self.resultsViewer.setObjectName("resultsViewer")
        self.resultsViewer.setColumnCount(0)
        self.resultsViewer.setRowCount(0)
        self.verticalLayout.addWidget(self.splitter)

        self.retranslateUi(SchemaSelectorWidget)
        QtCore.QMetaObject.connectSlotsByName(SchemaSelectorWidget)

    def retranslateUi(self, SchemaSelectorWidget):
        _translate = QtCore.QCoreApplication.translate
        SchemaSelectorWidget.setWindowTitle(_translate("SchemaSelectorWidget", "Form"))
