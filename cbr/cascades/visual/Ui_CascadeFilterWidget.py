# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'cbr/cascades/visual/CascadeFilterWidget.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_CascadeFilterWidget(object):
    def setupUi(self, CascadeFilterWidget):
        CascadeFilterWidget.setObjectName("CascadeFilterWidget")
        CascadeFilterWidget.resize(150, 42)
        self.horizontalLayout = QtWidgets.QHBoxLayout(CascadeFilterWidget)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.identitySpinBox = QtWidgets.QSpinBox(CascadeFilterWidget)
        self.identitySpinBox.setMaximum(100)
        self.identitySpinBox.setProperty("value", 80)
        self.identitySpinBox.setObjectName("identitySpinBox")
        self.horizontalLayout.addWidget(self.identitySpinBox)
        self.howComboBox = QtWidgets.QComboBox(CascadeFilterWidget)
        self.howComboBox.setObjectName("howComboBox")
        self.horizontalLayout.addWidget(self.howComboBox)

        self.retranslateUi(CascadeFilterWidget)
        QtCore.QMetaObject.connectSlotsByName(CascadeFilterWidget)

    def retranslateUi(self, CascadeFilterWidget):
        _translate = QtCore.QCoreApplication.translate
        CascadeFilterWidget.setWindowTitle(_translate("CascadeFilterWidget", "Form"))