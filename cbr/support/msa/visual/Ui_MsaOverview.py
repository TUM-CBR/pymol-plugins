# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/support/msa/visual/MsaOverview.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MsaOverview(object):
    def setupUi(self, MsaOverview):
        MsaOverview.setObjectName("MsaOverview")
        MsaOverview.resize(772, 478)
        self.verticalLayout = QtWidgets.QVBoxLayout(MsaOverview)
        self.verticalLayout.setObjectName("verticalLayout")
        self.msaOverviewScrollArea = QtWidgets.QScrollArea(MsaOverview)
        self.msaOverviewScrollArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.msaOverviewScrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.msaOverviewScrollArea.setWidgetResizable(True)
        self.msaOverviewScrollArea.setObjectName("msaOverviewScrollArea")
        self.msaOverviewWidget = QtWidgets.QWidget()
        self.msaOverviewWidget.setGeometry(QtCore.QRect(0, 0, 752, 458))
        self.msaOverviewWidget.setObjectName("msaOverviewWidget")
        self.msaOverviewScrollArea.setWidget(self.msaOverviewWidget)
        self.verticalLayout.addWidget(self.msaOverviewScrollArea)

        self.retranslateUi(MsaOverview)
        QtCore.QMetaObject.connectSlotsByName(MsaOverview)

    def retranslateUi(self, MsaOverview):
        _translate = QtCore.QCoreApplication.translate
        MsaOverview.setWindowTitle(_translate("MsaOverview", "Form"))